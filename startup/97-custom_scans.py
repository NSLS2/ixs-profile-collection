import bluesky.plans as bp
import bluesky.plan_stubs as bps
import bluesky.preprocessors as bpp
import numpy as np

# plt, np, db are IPython globals injected by nslsii.configure_base (00-startup.py).
# det1, det134, sclr, sr_curr are ophyd device globals from hardware startup files.
# set_lambda_exposure is defined in 96-lambda.py.
# _get_ah_scan_dets, _AH_DET_NAMES, short_label are defined in 97-setup_plans.py.
# select_detector_fields, CustomLiveMesh are imported via
#   "from utils.CustomLivePlot import *" in 00-startup.py.
# All of the above are available as IPython globals at the time any plan is called.

#*******************************************************************************************************
# Persistent figure for 2-D mesh scan display.
# plt and the Qt event loop are already initialised by 00-startup.py.
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
    Total events per scan: slow_steps × fast_steps.
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
    Total events per scan: slow_steps × fast_steps.
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
