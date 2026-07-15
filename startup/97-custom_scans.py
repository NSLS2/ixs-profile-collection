import bluesky.plans as bp
import bluesky.plan_stubs as bps
import bluesky.preprocessors as bpp
import numpy as np
import re
from bluesky.callbacks.fitting import PeakStats

# IPython globals available in all function bodies at call time:
#   plt, np, db          — injected by nslsii.configure_base (00-startup.py)
#   det1, det134, sclr, sr_curr
#                        — ophyd device globals from hardware startup files
#   set_lambda_exposure  — defined in 96-lambda.py
#   select_detector_fields, CustomLivePlot, CustomLiveMesh
#                        — star-imported from utils.CustomLivePlot in 00-startup.py

#*******************************************************************************************************
# Detector names eligible for ValuesPerRead/AveragingTime auto-configuration.
# Includes AH17x picoammeters (det1-det5) and TetrAMM (tm1).
_AH_DET_NAMES = frozenset({'det1', 'det2', 'det3', 'det4', 'det5', 'tm1'})


def _get_ah_scan_dets(primary_det):
    """Return the list of detector objects to configure for a scan.

    det1 is always included: det134 shares its hardware and participates
    in every ascan dets list; det1 is also always physically acquiring
    during dscan.

    The primary detector is appended if it is det2, det3, det4, det5,
    or tm1. If the primary detector is det1 it is already covered.
    tm2 is excluded and can be added to _AH_DET_NAMES if needed.
    """
    targets = [det1]
    if primary_det.name in (_AH_DET_NAMES - {'det1'}):
        targets.append(primary_det)
    return targets


#*******************************************************************************************************
# Persistent figure for 1-D scan display.
# plt and the Qt event loop are already initialised by 00-startup.py.
myfig, myaxs = plt.subplots(figsize=(8, 5), num="Live Scan", clear=False)
myfig.canvas.manager.set_window_title("Live Scan")
myfig.show()
myfig.canvas.draw_idle()
myfig.canvas.flush_events()

#*******************************************************************************************************
def short_label(field):
    """
    Convert long ophyd/bluesky field names into short legend mnemonics.

    Examples
    --------
    det2_current2_mean_value -> det2.2
    det3_current7_mean_value -> det3.7
    lambda_det_stats7_total  -> lambda.7
    lambda_det_stats3_total  -> lambda.3
    tm1_sum_all_mean_value    -> tm1.4   # special case
    tm2_sum_all_mean_value    -> tm2.4   # special case
    """
    if field in ("tm1_sum_all_mean_value", "tm2_sum_all_mean_value"):
        name = field.split("_", 1)[0]
        return f"{name}.4"

    # Generic current-channel detector pattern
    m = re.fullmatch(r"(det\d+)_current(\d+)_mean_value", field)
    if m:
        det, ch = m.groups()
        return f"{det}.{ch}"

    # tm1 / tm2 pattern
    m = re.fullmatch(r"(tm\d+)_current(\d+)_mean_value", field)
    if m:
        det_name, ch = m.groups()
        return f"{det_name}.{ch}"

    # Lambda stats pattern
    m = re.fullmatch(r"lambda_det_stats(\d+)_total", field)
    if m:
        ch = m.group(1)
        return f"lambda.{ch}"

    # Fallback: return original field name unchanged
    return field

#*******************************************************************************************************
def peaks_stats_print(dets_name, peak_stats):

#    headers = ["com","cen","fwhm","max","min"]
    print(dets_name)
#    print(peak_stats)
    COM = peak_stats['stats'].com
    if COM == None:
        print(f"COM: None")
    else:
        print(f"COM: {COM:.3f}")

    CEN = peak_stats['stats'].cen
    if CEN == None:
        print(f"CEN: None")
    else:
        print(f"CEN: {CEN:.3f}")

    FWHM = peak_stats['stats'].fwhm
    if FWHM == None:
        print(f"FWHM: None")
    else:
        print(f"FWHM: {FWHM:.3f}")
    pmax = peak_stats['stats'].max[1]
    if pmax < 1 or pmax > 1.e4:
        print(f"MAX: {pmax:.1e} at {peak_stats['stats'].max[0]:.3f}")
    else:
        print(f"MAX: {pmax:.1f} at {peak_stats['stats'].max[0]:.3f}")
    pmin = peak_stats['stats'].min[1]
    if pmin < 1 or pmin > 1.e4:
        print(f"MIN: {pmin:.1e} at {peak_stats['stats'].min[0]:.3f}")
    else:
        print(f"MIN: {pmin:.1f} at {peak_stats['stats'].min[0]:.3f}")


#*******************************************************************************************************
def dscan(mot, start, stop, steps, det, ct, det_ch=None, md=None):
    """
    Relative scan with live plotting and peak statistics.

    Parameters
    ----------
    mot : OphydObj
        Scanned motor.
    start, stop : float
        Relative scan limits.
    steps : int
        Number of points.
    det : OphydObj
        Main detector to plot.
    ct : float
        Count/exposure time.
    det_ch : list[int] or None
        Detector channels to plot.
        If None:
          - lambda_det -> default lambda_det_stats7_total
          - other detectors -> default channel 0
    md : dict or None
        Extra metadata.
    """
    md = md or {}
    md["count_time"] = ct

    if det_ch is None:
        det_ch = [0]

    dets = [det, det134, sclr.channels.chan13, sr_curr]
    # dets = [det, sclr.channels.chan13, sr_curr]

    # apply exposure if detector is lambda_det
    if getattr(det, "name", None) == "lambda_det":
        yield from set_lambda_exposure(ct)

    # --- AH501D picoammeter pre-configuration ---
    _ah_dets = _get_ah_scan_dets(det)
    _ah_old = {d.name: (d.values_per_read.get(), d.averaging_time.get())
               for d in _ah_dets}
    md["ah501d_config"] = {
        "values_per_read": 1,
        "averaging_time": ct,
        "detectors_configured": [d.name for d in _ah_dets],
    }
    for _d in _ah_dets:
        yield from bps.mv(_d.values_per_read, 1, _d.averaging_time, ct)

    # resolve actual event-data field names to plot/stat
    y_fields = select_detector_fields(det, det_ch)

    # compact legend labels
    legend_keys = {f: short_label(f) for f in y_fields}

    # one live plot callback for all selected fields
    liveplot_cb = CustomLivePlot(
        y_fields=y_fields,
        x=mot.name,
        ax=myaxs,
        legend_keys=legend_keys,
        clear_on_start=True,
        show_stats=True,
        update_every=1,
        title=f"{det.name} vs {mot.name}",
    )

    # one PeakStats per plotted field
    stats_list = [PeakStats(mot.name, field) for field in y_fields]

    subs_list = [liveplot_cb]
    subs_list.extend(stats_list)

    plan = bpp.subs_wrapper(
        bp.rel_scan(dets, mot, start, stop, steps, md=md),
        subs_list,
    )

    # --- AH501D restore wrapper (runs even on failure/interrupt) ---
    def _restore_ah():
        for _d in _ah_dets:
            _vpr, _avg = _ah_old[_d.name]
            yield from bps.mv(_d.values_per_read, _vpr, _d.averaging_time, _avg)

    plan = bpp.finalize_wrapper(plan, _restore_ah())

    yield from plan

    print("\n")
    for field, stats in zip(y_fields, stats_list):
        peaks_stats_print(field, stats)
        print("\n")

    return stats_list

#*******************************************************************************************************
def ascan(mot, start, stop, steps, det, ct, det_ch=None, md=None):
    """
    Absolute scan with live plotting and peak statistics.

    Parameters
    ----------
    mot : OphydObj
        Scanned motor.
    start, stop : float
        Absolute scan limits.
    steps : int
        Number of points.
    det : OphydObj
        Main detector to plot.
    ct : float
        Count/exposure time.
    det_ch : int, list[int], or None
        Detector channels to plot.
        If None:
          - lambda_det -> default lambda_det_stats7_total
          - other detectors -> default channel 0
    md : dict or None
        Extra metadata.
    """
    md = md or {}
    md["count_time"] = ct

    if det_ch is None:
        det_ch = [0]

    dets = [det, det134, sclr.channels.chan13, sr_curr]
    # dets = [det]

    # apply exposure if detector is lambda_det
    if getattr(det, "name", None) == "lambda_det":
        yield from set_lambda_exposure(ct)

    # --- AH501D picoammeter pre-configuration ---
    _ah_dets = _get_ah_scan_dets(det)
    _ah_old = {d.name: (d.values_per_read.get(), d.averaging_time.get())
               for d in _ah_dets}
    md["ah501d_config"] = {
        "values_per_read": 1,
        "averaging_time": ct,
        "detectors_configured": [d.name for d in _ah_dets],
    }
    for _d in _ah_dets:
        yield from bps.mv(_d.values_per_read, 1, _d.averaging_time, ct)

    # resolve actual event-data field names to plot/stat
    y_fields = select_detector_fields(det, det_ch)
    if not y_fields:
        raise RuntimeError(f"No plot fields resolved for detector {det.name}")

    # compact legend labels
    legend_keys = {f: short_label(f) for f in y_fields}

    # one live plot callback for all selected fields
    liveplot_cb = CustomLivePlot(
        y_fields=y_fields,
        x=mot.name,
        ax=myaxs,
        legend_keys=legend_keys,
        clear_on_start=True,
        show_stats=True,
        update_every=1,
        title=f"{det.name} vs {mot.name}",
    )

    # one PeakStats per plotted field
    stats_list = [PeakStats(mot.name, field) for field in y_fields]

    subs_list = [liveplot_cb]
    subs_list.extend(stats_list)

    plan = bpp.subs_wrapper(
        bp.scan(dets, mot, start, stop, steps, md=md),
        subs_list,
    )

    # --- AH501D restore wrapper (runs even on failure/interrupt) ---
    def _restore_ah():
        for _d in _ah_dets:
            _vpr, _avg = _ah_old[_d.name]
            yield from bps.mv(_d.values_per_read, _vpr, _d.averaging_time, _avg)

    plan = bpp.finalize_wrapper(plan, _restore_ah())

    yield from plan

    print("\n")
    for field, stats in zip(y_fields, stats_list):
        peaks_stats_print(field, stats)
        print("\n")

    return stats_list

#*******************************************************************************************************
# Persistent figure for 2-D mesh scan display.
mymeshfig = plt.figure(figsize=(10, 6), num="Live Mesh Scan", clear=False)
mymeshfig.canvas.manager.set_window_title("Live Mesh Scan")
mymeshfig.show()
mymeshfig.canvas.draw_idle()
mymeshfig.canvas.flush_events()

#*******************************************************************************************************
def mesh_stats_print(field, slow_name, fast_name, table):
    """
    Print the maximum value and its (slow, fast) motor position for one mesh scan field.

    Parameters
    ----------
    field : str
        Event-data field name (column in *table*).
    slow_name : str
        Column name for the slow (outer) motor position in *table*.
    fast_name : str
        Column name for the fast (inner) motor position in *table*.
    table : pandas.DataFrame
        Table returned by db[-1].table().
    """
    print(field)
    if field not in table.columns:
        print("  (field not found in scan table)")
        return

    idx = table[field].idxmax()
    value = table[field][idx]
    slow_val = table[slow_name][idx] if slow_name in table.columns else float("nan")
    fast_val = table[fast_name][idx] if fast_name in table.columns else float("nan")

    if value < 1 or value > 1.0e4:
        print(f"MAX: {value:.4e} at {slow_name}={slow_val:.5g}, {fast_name}={fast_val:.5g}")
    else:
        print(f"MAX: {value:.4f} at {slow_name}={slow_val:.5g}, {fast_name}={fast_val:.5g}")


#*******************************************************************************************************
def dmeshscan(slow_mot, slow_start, slow_stop, slow_steps,
              fast_mot, fast_start, fast_stop, fast_steps,
              det, ct, det_ch=None, md=None):
    """
    Relative 2-D mesh scan (bp.rel_grid_scan) with live scatter plot and
    post-scan statistics.

    Parameters
    ----------
    slow_mot : ophyd device
        Outer (slow) motor — steps once per complete inner sweep.
    slow_start, slow_stop : float
        Relative scan limits for the slow axis (offset from current position).
    slow_steps : int
        Number of points on the slow axis.
    fast_mot : ophyd device
        Inner (fast) motor — sweeps fully at every slow-axis position.
    fast_start, fast_stop : float
        Relative scan limits for the fast axis.
    fast_steps : int
        Number of points on the fast axis.
    det : ophyd device
        Primary detector.
    ct : float
        Count / exposure time in seconds.
    det_ch : list[int] or None
        Indices into det.hints['fields'] to plot.  Defaults to [0].
    md : dict or None
        Extra metadata merged into the run start document.

    Returns
    -------
    None
        Post-scan statistics are printed to stdout.

    Notes
    -----
    Total events per scan: slow_steps x fast_steps.
    snake_axes is hardcoded False (_SNAKE_AXES constant).  Change it to True
    to enable boustrophedon scanning without altering the function signature.
    """
    md = md or {}
    md["count_time"] = ct

    if det_ch is None:
        det_ch = [0]

    dets = [det, det134, sclr.channels.chan13, sr_curr]

    # Lambda special case: programme exposure time before AH501D config
    if getattr(det, "name", None) == "lambda_det":
        yield from set_lambda_exposure(ct)

    # --- AH501D picoammeter pre-configuration ---
    _ah_dets = _get_ah_scan_dets(det)
    _ah_old = {d.name: (d.values_per_read.get(), d.averaging_time.get())
               for d in _ah_dets}
    md["ah501d_config"] = {
        "values_per_read": 1,
        "averaging_time": ct,
        "detectors_configured": [d.name for d in _ah_dets],
    }
    for _d in _ah_dets:
        yield from bps.mv(_d.values_per_read, 1, _d.averaging_time, ct)

    # Field and subplot-title resolution
    y_fields = select_detector_fields(det, det_ch)
    legend_keys = {f: short_label(f) for f in y_fields}

    # Live 2-D scatter callback
    meshplot_cb = CustomLiveMesh(
        y_fields=y_fields,
        outer_field=slow_mot.name,
        inner_field=fast_mot.name,
        fig=mymeshfig,
        legend_keys=legend_keys,
        clear_on_start=True,
        update_every=10,
        title=f"{det.name}: {slow_mot.name} vs {fast_mot.name}",
    )

    _SNAKE_AXES = False
    plan = bpp.subs_wrapper(
        bp.rel_grid_scan(
            dets,
            slow_mot, slow_start, slow_stop, slow_steps,
            fast_mot, fast_start, fast_stop, fast_steps,
            _SNAKE_AXES,
            md=md,
        ),
        [meshplot_cb],
    )

    # --- AH501D restore wrapper (runs on completion, failure, or Ctrl+C) ---
    def _restore_ah():
        for _d in _ah_dets:
            _vpr, _avg = _ah_old[_d.name]
            yield from bps.mv(_d.values_per_read, _vpr, _d.averaging_time, _avg)

    plan = bpp.finalize_wrapper(plan, _restore_ah())
    yield from plan

    # Post-scan statistics
    print("\n")
    table = db[-1].table()
    for field in y_fields:
        mesh_stats_print(field, slow_mot.name, fast_mot.name, table)
        print("\n")


#*******************************************************************************************************
def ameshscan(slow_mot, slow_start, slow_stop, slow_steps,
              fast_mot, fast_start, fast_stop, fast_steps,
              det, ct, det_ch=None, md=None):
    """
    Absolute 2-D mesh scan (bp.grid_scan) with live scatter plot and
    post-scan statistics.

    Parameters
    ----------
    slow_mot : ophyd device
        Outer (slow) motor — steps once per complete inner sweep.
    slow_start, slow_stop : float
        Absolute positions for the slow axis.
    slow_steps : int
        Number of points on the slow axis.
    fast_mot : ophyd device
        Inner (fast) motor — sweeps fully at every slow-axis position.
    fast_start, fast_stop : float
        Absolute positions for the fast axis.
    fast_steps : int
        Number of points on the fast axis.
    det : ophyd device
        Primary detector.
    ct : float
        Count / exposure time in seconds.
    det_ch : list[int] or None
        Indices into det.hints['fields'] to plot.  Defaults to [0].
    md : dict or None
        Extra metadata merged into the run start document.

    Returns
    -------
    None
        Post-scan statistics are printed to stdout.

    Notes
    -----
    Total events per scan: slow_steps x fast_steps.
    snake_axes is hardcoded False (_SNAKE_AXES constant).
    Raises RuntimeError if no fields resolve for the given detector and det_ch.
    """
    md = md or {}
    md["count_time"] = ct

    if det_ch is None:
        det_ch = [0]

    dets = [det, det134, sclr.channels.chan13, sr_curr]

    # Lambda special case: programme exposure time before AH501D config
    if getattr(det, "name", None) == "lambda_det":
        yield from set_lambda_exposure(ct)

    # --- AH501D picoammeter pre-configuration ---
    _ah_dets = _get_ah_scan_dets(det)
    _ah_old = {d.name: (d.values_per_read.get(), d.averaging_time.get())
               for d in _ah_dets}
    md["ah501d_config"] = {
        "values_per_read": 1,
        "averaging_time": ct,
        "detectors_configured": [d.name for d in _ah_dets],
    }
    for _d in _ah_dets:
        yield from bps.mv(_d.values_per_read, 1, _d.averaging_time, ct)

    # Field and subplot-title resolution
    y_fields = select_detector_fields(det, det_ch)
    if not y_fields:
        raise RuntimeError(f"No plot fields resolved for detector {det.name!r}")

    legend_keys = {f: short_label(f) for f in y_fields}

    # Live 2-D scatter callback
    meshplot_cb = CustomLiveMesh(
        y_fields=y_fields,
        outer_field=slow_mot.name,
        inner_field=fast_mot.name,
        fig=mymeshfig,
        legend_keys=legend_keys,
        clear_on_start=True,
        update_every=10,
        title=f"{det.name}: {slow_mot.name} vs {fast_mot.name}",
    )

    _SNAKE_AXES = False
    plan = bpp.subs_wrapper(
        bp.grid_scan(
            dets,
            slow_mot, slow_start, slow_stop, slow_steps,
            fast_mot, fast_start, fast_stop, fast_steps,
            _SNAKE_AXES,
            md=md,
        ),
        [meshplot_cb],
    )

    # --- AH501D restore wrapper (runs on completion, failure, or Ctrl+C) ---
    def _restore_ah():
        for _d in _ah_dets:
            _vpr, _avg = _ah_old[_d.name]
            yield from bps.mv(_d.values_per_read, _vpr, _d.averaging_time, _avg)

    plan = bpp.finalize_wrapper(plan, _restore_ah())
    yield from plan

    # Post-scan statistics
    print("\n")
    table = db[-1].table()
    for field in y_fields:
        mesh_stats_print(field, slow_mot.name, fast_mot.name, table)
        print("\n")
