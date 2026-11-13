import bluesky.plans as bp
import bluesky.plan_stubs as bps
import bluesky.preprocessors as bpp
from bluesky.utils import RequestAbort, RequestStop
import numpy as np
import pandas as pd
import lmfit
from bluesky.callbacks import LiveFit
from bluesky.callbacks.mpl_plotting import LiveGrid
from tabulate import tabulate
from utils.mcm_alignment import (
    McmAlignConfig, CorrectionLimitError, HRM_IN, HRM_OUT, CONVERGED, CORRECT,
    classify_hrm_state, hrm_measurement_settings, decide_analyzer_move,
    analyzer_readback_ok, tx_correction, dch_moves, dy_correction,
    convergence_step)

# _AH_DET_NAMES, _get_ah_scan_dets, myfig, short_label, peaks_stats_print,
# dscan, ascan are defined in 97-custom_scans.py (loads before this file).

#*******************************************************************************************************
def gaussian(x, A, sigma, x0):
    exponent = -(x - x0)**2/(2 * sigma**2)
    # Clip exponent to prevent overflow; exp(-745) ≈ 0, exp(709) ≈ max float
    exponent = np.clip(exponent, -745, 0)
    return A*np.exp(exponent)


#*******************************************************************************************************
def stepup(x, A, sigma, x0, b):
    return A*(1-1/(1+np.exp((x-x0)/sigma)))+b


#*******************************************************************************************************
def stepdown(x, A, sigma, x0, b):
    return A*(1-1/(1+np.exp(-(x-x0)/sigma)))+b


#*******************************************************************************************************
def calc_lmfit(uid=-1, x="hrmE", channel=7):
    # Calculates fitting parameters for Gaussian function for energy scan with UID and Lambda channel
    y = f'lambda_det_stats{channel}_total'
    run = trc[uid]
    table = run['primary'].read([x, y])
    model = lmfit.Model(gaussian)
    
    lf = LiveFit(model, y, {'x': x}, {'A': table[y].max(), 'sigma': 0.7, 'x0': table[x][np.argmax(table[y].data)+1]})
    for name, doc in run.documents():
        lf(name, doc)
    gauss = gaussian(table[x], **lf.result.values)
    _ax = myfig.gca()
    _ax.plot(table[x], table[y], label=f"raw, channel={channel}", marker = 'o', linestyle = 'none')
    _ax.plot(table[x], gauss.values, label=f"gaussian fit {channel}")
    _ax.legend()

    myfig.canvas.draw_idle()
    myfig.canvas.flush_events()

    return lf.result.values


#*******************************************************************************************************
def calc_stepup_fit(x):
    # Calculates fitting parameters for step up function for MCM slits scan
    y = 'det2_current1_mean_value'
    run = trc[-1]
    table = run['primary'].read([x, y])
    model = lmfit.Model(stepup)
    lf = LiveFit(model, y, {'x': x}, {'A': table[y].max(), 'sigma': 0.25, 'x0': 0, 'b':0})
    for name, doc in run.documents():
        lf(name, doc)
    print(lf.result.values)
    stup = stepup(table[x], **lf.result.values)
    plt.clf()
    plt.plot(table[x], table[y], label=f"raw data", marker = 'o', linestyle = 'none')
    plt.plot(table[x], stup.values, label=f"data fit")
    plt.legend()
    return lf.result.values['x0']


#*******************************************************************************************************
def calc_stepdwn_fit(x):
    # Calculates fitting parameters for step down function for MCM slits scan
    y = 'det2_current1_mean_value'
    run = trc[-1]
    table = run['primary'].read([x, y])
    model = lmfit.Model(stepdown)
    lf = LiveFit(model, y, {'x': x}, {'A': table[y].max(), 'sigma': 0.25, 'x0': 0, 'b':0})
    for name, doc in run.documents():
        lf(name, doc)
    print(lf.result.values)
    stdw = stepdown(table[x], **lf.result.values)
    plt.clf()
    plt.plot(table[x], table[y], label=f"raw data", marker = 'o', linestyle = 'none')
    plt.plot(table[x], stdw.values, label=f"data fit")
    plt.legend()
    return lf.result.values['x0']


#*******************************************************************************************************
def GCarbon_Qscan(exp_time=2):
    # Test plan for the energy resolution at Q=1.2 with the Glassy Carbon
    Qq = [1.2]
    yield from bps.mv(analyzer_slits.top, 1, analyzer_slits.bottom, -1, analyzer_slits.outboard, 1.5, analyzer_slits.inboard, -1.5)
    yield from bps.mv(anapd, 25, whl, 0)
#    myplt = plotselect('lambda_det_stats7_total', hrmE.name)

    for kk in range(1):
        for q in Qq:
            th = qq2th(q)
            yield from bps.mv(spec.tth, th)
            # yield from set_lambda_exposure(exp_time)
            yield from dscan(hrmE, -10, 10, 100, lambda_det, exp_time)


#*******************************************************************************************************
def DxtalTempCalc(uid=-1):
    """
    Calculates temperature correction for the D crystals.
    
    Parameters
    ----------
    uid : int, optional
          relative scan position in the database
    
    """
    
    E0 = 9131.7     # energy (eV)
    TH = 88.5       # Dxtal asymmetry angle (deg)
    C1 = 3.725e-6   # constant (1/K)
    C2 = 5.88e-3    # constant (1/K)
    C3 = 5.548e-10  # constant (1/K2)
    T1 = 124.0      # temperature (K)
    T0 = 300.15     # crystal average temperature (K)

    bet = C1*(1 - np.exp(-C2*(T0-T1))) + C3*T0
    dE = []
    for n in range(1,7):
        fit_par = calc_lmfit(uid, channel=n)
        if fit_par['A'] < 100:
            print('**********************************')
            print('         WARNING !')
            print('      Fitting Error')
            return
        
        dE.append(fit_par['x0'])
    
    dE = [x-dE[0] for x in dE]
    dTe = [1.e-3*x/E0/bet for x in dE]
    dTh = [-1.e3*x*np.tan(np.radians(TH))/E0 for x in dE]
    
    DTe = [ura_temp.d1temp.read()['uratemperature_d1temp']['value']+dTe[0], 
           ura_temp.d2temp.read()['uratemperature_d2temp']['value']+dTe[1], 
           ura_temp.d3temp.read()['uratemperature_d3temp']['value']+dTe[2], 
           ura_temp.d4temp.read()['uratemperature_d4temp']['value']+dTe[3], 
           ura_temp.d5temp.read()['uratemperature_d5temp']['value']+dTe[4], 
           ura_temp.d6temp.read()['uratemperature_d6temp']['value']+dTe[5]]
    Dheader = [' ', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6']
    dE.insert(0,'dEnrg')
    dTe.insert(0,'dTemp')
    dTh.insert(0,'dThe')
    DTe.insert(0,'Dtemp')
    Ddata = [dE, dTh, dTe, DTe]
    print('---------------------------------------------------------------------')
    print(tabulate(Ddata, headers=Dheader, tablefmt='pipe', stralign='center', floatfmt='.4f'))
    print('---------------------------------------------------------------------\n')
    print('Select from the following options:')
    print('1. Update temperatures')
    print('2. Update Dxtals angles')
    print('3. Exit without updates')
    update_opts = input('Your choice: ')
    if update_opts == '1':
        ura_temp.d1temp.set(DTe[1])
        ura_temp.d2temp.set(DTe[2])
        ura_temp.d3temp.set(DTe[3])
        ura_temp.d4temp.set(DTe[4])
        ura_temp.d5temp.set(DTe[5])
        ura_temp.d6temp.set(DTe[6])
        print('\n')
        print('The temperatures are updated\n')
    elif update_opts == '2':
        yield from bps.mvr(analyzer_xtals.d2the, dTh[2], analyzer_xtals.d3the, dTh[3], analyzer_xtals.d4the, dTh[4], analyzer_xtals.d5the, dTh[5], analyzer_xtals.d6the, dTh[6],)
        print('\n')
        print('The Dxtals angles are updated\n')
    else:
        print('\n')
        print('Update is canceled\n')
#    return {'dEn':dE, 'dTem':dTe, 'dThe':dTh, 'DTem':DTe}


#*******************************************************************************************************
def ugap_setup():
#   Scans the ID gap and sets it to max
#
    old_gap = ivu22.position
    res = yield from dscan(ivu22, -20, 20, 21, tm1, 1, det_ch=[4])
    
    max_pos = res[0].max
    x_pos = max_pos[0]
    yield from bps.mv(ivu22, x_pos)
    print('\n')
    print(f'ID gap moved from {old_gap:.2f} um to {x_pos:.2f} um')
    print('ID gap alignment finished\n')


#*******************************************************************************************************
def _crl_validate(stats, initial_pos):
    """Return (xmax, contrast); raise AlignmentFailure with diagnostics."""
    MAX_SHIFT = 0.1
    MIN_CONTRAST = 0.2
    MIN_INTENSITY = 3.e-5
    fmax, fmin = stats.max, stats.min
    contrast = None
    try:
        ymax, ymin = fmax[1], fmin[1]
        if (ymax + ymin) != 0:
            contrast = (ymax - ymin) / (ymax + ymin)
    except (TypeError, IndexError):
        pass
    diag = f'initial={initial_pos:.3f} max={fmax} contrast={contrast}'
    if fmax is None or contrast is None:
        raise AlignmentFailure(f'CRL y scan found no intensity maximum ({diag}).')
    xmax = fmax[0]
    if not fmax[1] > MIN_INTENSITY:
        raise AlignmentFailure(
            f'CRL y maximum intensity {fmax[1]:.3e} not above {MIN_INTENSITY:.1e} ({diag}).')
    if abs(xmax - initial_pos) > MAX_SHIFT:
        raise AlignmentFailure(
            f'CRL y maximum {xmax - initial_pos:+.3f} mm from initial exceeds '
            f'+/-{MAX_SHIFT} mm ({diag}).')
    if contrast < MIN_CONTRAST:
        raise AlignmentFailure(
            f'CRL y contrast {contrast:.2f} below {MIN_CONTRAST} ({diag}).')
    return xmax, contrast


def crl_setup():
    """Move crl.y to the intensity maximum. Raises AlignmentFailure on failure.

    Failure restores crl.y to its initial position and aborts the calling plan.
    """
    initial_pos = crl.y.position
    status, reason, restore_errors = 'success', '', []
    xmax = contrast = None
    try:
        res = yield from dscan(crl.y, -0.2, 0.2, 20, tm1, 1, det_ch=[4])
        xmax, contrast = _crl_validate(res[0], initial_pos)
        yield from bps.mv(crl.y, xmax)
    except (RequestAbort, RequestStop):
        raise
    except Exception as exc:
        status, reason = 'failed', _mcm_describe(exc)
        try:
            yield from bps.mv(crl.y, initial_pos)
        except (RequestAbort, RequestStop):
            raise
        except Exception as rexc:
            restore_errors.append(f'crl.y: {_mcm_describe(rexc)}')
    try:
        final = crl.y.position
    except Exception as exc:
        final = f'unreadable ({exc!r})'
    print('\n' + '*' * 60)
    print(f'CRL setup: {status.upper()}' + (f' - {reason}' if reason else ''))
    print(f'  initial: {initial_pos}  final: {final}')
    print(f'  max position: {xmax}  contrast: {contrast}')
    print(f'  restore_errors: {restore_errors}')
    print('*' * 60 + '\n')
    if status == 'failed':
        raise AlignmentFailure(reason)


#*******************************************************************************************************
def _dcm_validate(stats, initial_pos):
    """Return (cen, fwhm, contrast); raise AlignmentFailure with diagnostics."""
    FWHM_MIN, FWHM_MAX = 46.0, 47.0
    CONTRAST_REF, CONTRAST_TOL = 0.9, 0.1
    MAX_SHIFT = 20.0
    MIN_INTENSITY = 3.e-5
    fwhm, cen = stats.fwhm, stats.cen
    contrast = None
    try:
        ymax, ymin = stats.max[1], stats.min[1]
        if (ymax + ymin) != 0:
            contrast = (ymax - ymin) / (ymax + ymin)
    except (TypeError, IndexError):
        pass
    diag = f'initial={initial_pos:.3f} fwhm={fwhm} cen={cen} contrast={contrast}'
    if fwhm is None or cen is None or contrast is None:
        raise AlignmentFailure(f'DCM p1 scan returned invalid results ({diag}).')
    if not stats.max[1] > MIN_INTENSITY:
        raise AlignmentFailure(
            f'DCM p1 maximum intensity {stats.max[1]:.3e} not above '
            f'{MIN_INTENSITY:.1e} ({diag}).')
    if not FWHM_MIN <= fwhm <= FWHM_MAX:
        raise AlignmentFailure(
            f'DCM p1 FWHM {fwhm:.2f} outside [{FWHM_MIN}, {FWHM_MAX}] ({diag}).')
    if abs(contrast - CONTRAST_REF) > CONTRAST_TOL + 1e-9:
        raise AlignmentFailure(
            f'DCM p1 contrast {contrast:.2f} outside '
            f'{CONTRAST_REF}+/-{CONTRAST_TOL} ({diag}).')
    if abs(cen - initial_pos) > MAX_SHIFT:
        raise AlignmentFailure(
            f'DCM p1 shift {cen - initial_pos:+.2f} exceeds +/-{MAX_SHIFT} ({diag}).')
    return cen, fwhm, contrast


def dcm_setup():
    """Center dcm.p1 on the beam-intensity peak. Raises AlignmentFailure on failure.

    Failure restores dcm.p1 to its initial position and aborts the calling plan.
    Operator interruption (RequestAbort/RequestStop) restores nothing.
    """
    initial_pos = dcm.p1.position
    status, reason, restore_errors = 'success', '', []
    cen = fwhm = contrast = None
    try:
        res = yield from dscan(dcm.p1, -80, 80, 40, tm1, 1, det_ch=[4])
        cen, fwhm, contrast = _dcm_validate(res[0], initial_pos)
        yield from bps.mv(dcm.p1, cen)
    except (RequestAbort, RequestStop):
        raise
    except OperatorAbort as exc:
        status, reason = 'aborted', str(exc)
    except Exception as exc:
        status, reason = 'failed', _mcm_describe(exc)
        try:
            yield from bps.mv(dcm.p1, initial_pos)
        except (RequestAbort, RequestStop):
            raise
        except Exception as rexc:
            restore_errors.append(f'dcm.p1: {_mcm_describe(rexc)}')
    try:
        final = dcm.p1.position
    except Exception as exc:
        final = f'unreadable ({exc!r})'
    print('\n' + '*' * 60)
    print(f'DCM setup: {status.upper()}' + (f' - {reason}' if reason else ''))
    print(f'  initial: {initial_pos}  final: {final}')
    print(f'  fwhm: {fwhm}  contrast: {contrast}  cen: {cen}')
    print(f'  restore_errors: {restore_errors}')
    print('*' * 60 + '\n')
    if status == 'failed':
        raise AlignmentFailure(reason)
    if status == 'aborted':
        raise OperatorAbort(reason)


#*******************************************************************************************************
# MCM alignment. Pure logic (HRM state, corrections, limits, convergence) is in
# utils/mcm_alignment.py. Requirements: docs/mcm_alignment_scope.md.
class AlignmentFailure(RuntimeError):
    """Alignment failed. Triggers failure restoration."""


class OperatorAbort(RuntimeError):
    """Operator declined a prompt. No restoration is performed."""


def _mcm_describe(exc):
    """Readable failure text; includes the cause of a failed device status."""
    text = str(exc) or repr(exc)
    try:
        inner = exc.args[0].exception()
    except Exception:
        inner = None
    return f'{text}: {inner!r}' if inner is not None else text


def _mcm_yes(answer):
    return str(answer).strip().lower() in ('y', 'yes')


def _mcm_travel(motor):
    """Motor (low, high) travel limits, or None if unavailable or unset."""
    try:
        low, high = motor.limits
    except Exception:
        return None
    return None if low == high else (low, high)


def _mcm_new_state(cfg):
    return {
        'cfg': cfg, 'hrm_state': None, 'hrm_readbacks': None,
        'saved': {}, 'touched': set(), 'prep_started': False,
        'mcm_ypos': None, 'mcm_x_active': False, 'mcm_y_return': None,
        'vertical_started': False,
        'axes': {}, 'history': [], 'corrections': [], 'uids': [],
        'restored': False, 'restore_errors': [],
    }


def mcm_setup_prep(cfg, state):
    """Save entry values and prepare the measurement (scope section 5)."""
    try:
        ux, dx = hrm2.ux.position, hrm2.dx.position
    except Exception:
        ux = dx = None
    hrm_state = classify_hrm_state(ux, dx, cfg)
    state['hrm_state'], state['hrm_readbacks'] = hrm_state, (ux, dx)
    if hrm_state not in (HRM_IN, HRM_OUT):
        raise AlignmentFailure(
            f'HRM state is {hrm_state} (hrm2.ux = {ux}, hrm2.dx = {dx}; in = '
            f'{cfg.hrm_in_position:g} +/- {cfg.hrm_tolerance:g} mm, out = '
            f'{cfg.hrm_out_position:g} +/- {cfg.hrm_tolerance:g} mm). '
            'Nothing was moved.')

    saved = state['saved']
    saved['y'] = anc_xtal.y.position
    saved['slits'] = {
        'top': analyzer_slits.top.position,
        'bottom': analyzer_slits.bottom.position,
        'outboard': analyzer_slits.outboard.position,
        'inboard': analyzer_slits.inboard.position,
    }
    saved['anpd'] = anpd.position
    saved['whl'] = whl.position
    saved['em_range'] = det2.em_range.get()
    if hrm_state == HRM_IN:
        state['mcm_ypos'] = mcm.y.position  # MCM_YPOS, fixed for the run

    decision = decide_analyzer_move(saved['y'], cfg)
    print(decision.message)
    if not decision.ok:
        raise AlignmentFailure(decision.message)

    state['prep_started'] = True
    if decision.move:
        state['touched'].add('y')
        yield from bps.mv(anc_xtal.y, decision.y_target)
        readback = anc_xtal.y.position
        if not analyzer_readback_ok(readback, cfg):
            raise AlignmentFailure(
                f'Analyzer readback {readback} mm is outside '
                f'({cfg.analyzer_min:g}, {cfg.analyzer_max:g}) mm after moving '
                f'from {decision.y_initial} to {decision.y_target} mm.')

    state['touched'].add('slits')
    yield from bps.mv(analyzer_slits.top, cfg.slit_top,
                      analyzer_slits.bottom, cfg.slit_bottom,
                      analyzer_slits.outboard, cfg.slit_outboard,
                      analyzer_slits.inboard, cfg.slit_inboard)
    state['touched'].add('anpd')
    yield from bps.mv(anpd, cfg.anpd_position)
    settings = hrm_measurement_settings(hrm_state, cfg)
    state['touched'].update(('whl', 'em_range'))
    yield from bps.mv(whl, settings['whl'], det2.em_range, settings['em_range'])
    # airpad is not restored; spec.tth is never restored.
    yield from bps.mv(airpad, cfg.airpad_on)
    yield from bps.mv(spec.tth, cfg.spec_tth_position)

    intensity = yield from bps.rd(det2.current1.mean_value)
    state['prep_intensity'] = intensity
    if intensity < cfg.low_intensity_threshold:
        raise AlignmentFailure(
            f'Low intensity on det2.current1: {intensity:.3g} < '
            f'{cfg.low_intensity_threshold:.3g}.')


def mcm_setup_post(state):
    """Restore the saved entry values that preparation changed.

    Every item is attempted. Errors are recorded in state['restore_errors']
    and never raised, so they cannot hide an original failure.
    """
    saved, touched = state['saved'], state['touched']
    targets = []
    if 'y' in touched:
        targets.append(('anc_xtal.y', (anc_xtal.y, saved['y'])))
    if 'slits' in touched:
        s = saved['slits']
        targets.append(('analyzer_slits', (
            analyzer_slits.top, s['top'], analyzer_slits.bottom, s['bottom'],
            analyzer_slits.outboard, s['outboard'],
            analyzer_slits.inboard, s['inboard'])))
    if 'anpd' in touched:
        targets.append(('anpd', (anpd, saved['anpd'])))
    if 'whl' in touched:
        targets.append(('whl', (whl, saved['whl'])))
    if 'em_range' in touched:
        targets.append(('det2.em_range', (det2.em_range, saved['em_range'])))
    for label, args in targets:
        try:
            yield from bps.mv(*args)
        except (RequestAbort, RequestStop):
            raise
        except Exception as exc:
            state['restore_errors'].append(f'{label}: {_mcm_describe(exc)}')
    state['restored'] = not state['restore_errors']


def _mcm_scan_fit(motor, cfg, state):
    """Scan `motor` (relative, centered on its current position) and fit it.

    Uses the exact scan UID, flushing the Tiled writer after the scan and
    before the run is read back. Returns (position of maximum, fitted
    intensity, uid).
    """
    uid_holder = {}

    def _capture(name, doc):
        uid_holder['uid'] = doc['uid']

    plan = dscan(motor, -cfg.scan_half_range, cfg.scan_half_range,
                 cfg.scan_points, det2, cfg.scan_time, det_ch=[0])
    yield from bpp.subs_wrapper(plan, {'start': [_capture]})
    uid = uid_holder.get('uid')
    if uid is None:
        raise AlignmentFailure(f'No Start document captured for the {motor.name} scan.')
    state['uids'].append(uid)
    try:
        flush_tiled()
        pos, intensity = calculate_max_value(
            uid=uid, x=motor.name, y='det2_current1_mean_value',
            delta=cfg.fit_delta, sampling=cfg.fit_sampling)
    except Exception as exc:
        raise AlignmentFailure(
            f'Fit failed for {motor.name} scan {uid}: {exc}') from exc
    return float(pos), float(intensity), uid


def _mcm_confirm(ask, cfg, motor, deltas):
    """Ask before a correction (interactive mode). Declining aborts."""
    if not cfg.interactive:
        return
    pos = motor.position
    target = pos + sum(deltas)
    steps = ', '.join(f'{d:+.6g}' for d in deltas)
    print(f'Proposed: {motor.name} at {pos:.6g}, relative move(s) {steps}, target {target:.6g}')
    if not _mcm_yes(ask('Apply? (yes/no): ')):
        raise OperatorAbort(f'Correction of {motor.name} declined by operator.')


def _mcm_correction_moves(axis, hrm_state, err, cfg):
    """Return (motor, [relative moves]) after limit checks; never clips."""
    try:
        if axis == 'horizontal':
            if hrm_state == HRM_OUT:
                m = sample_stage.tx
                return m, [tx_correction(err, cfg, m.position, _mcm_travel(m))]
            m = hrm2.dch
            return m, dch_moves(err, cfg, m.position, _mcm_travel(m))
        m = hrm2.dy
        return m, [dy_correction(err, cfg, m.position, _mcm_travel(m))]
    except CorrectionLimitError as exc:
        raise AlignmentFailure(str(exc)) from exc


def _mcm_axis_loop(axis, motor, reference, tolerance, cfg, state, ask):
    """Scan, compare with `reference`, correct and rescan (max_scans)."""
    for n in range(1, cfg.max_scans + 1):
        yield from bps.mv(motor, reference)
        pos, intensity, uid = yield from _mcm_scan_fit(motor, cfg, state)
        yield from bps.mv(motor, reference)
        err = pos - reference
        state['history'].append({
            'axis': axis, 'scan': n, 'uid': uid, 'max': pos,
            'reference': reference, 'error': err, 'intensity': intensity})
        decision = convergence_step(n, err, tolerance, cfg.max_scans)
        print(f'{axis} scan {n}: max {pos:.5f}, reference {reference:.5f}, error {err:+.5f}')
        if decision == CONVERGED:
            return
        if decision != CORRECT:
            raise AlignmentFailure(
                f'{axis} alignment did not converge after {n} scans '
                f'(error {err:+.5f}, tolerance {tolerance:g}).')
        corr_motor, moves = _mcm_correction_moves(axis, state['hrm_state'], err, cfg)
        _mcm_confirm(ask, cfg, corr_motor, moves)
        for delta in moves:
            yield from bps.mvr(corr_motor, delta)
        state['corrections'].append({
            'axis': axis, 'scan': n, 'motor': corr_motor.name, 'moves': list(moves)})


def _mcm_axis_enabled(axis, enabled, cfg, ask):
    if not enabled:
        return False
    if cfg.interactive:
        return _mcm_yes(ask(f'Perform {axis} alignment? (yes/no): '))
    return True


def _mcm_vertical_hrm_out(cfg, state, ask):
    pre = mcm.y.position
    state['mcm_y_return'] = pre
    ymax, intensity, uid = yield from _mcm_scan_fit(mcm.y, cfg, state)
    state['history'].append({
        'axis': 'vertical', 'scan': 1, 'uid': uid, 'max': ymax,
        'reference': pre, 'error': ymax - pre, 'intensity': intensity})
    print(f'vertical scan: max {ymax:.5f}, pre-scan position {pre:.5f}')
    _mcm_confirm(ask, cfg, mcm.y, [ymax - pre])
    yield from bps.mv(mcm.y, ymax)
    state['mcm_y_return'] = None  # moved on purpose; nothing to return to


def _mcm_failure_cleanup(state):
    """Return MCM axes and restore preparation after an alignment failure."""
    cfg = state['cfg']

    def _try(label, plan):
        try:
            yield from plan
        except (RequestAbort, RequestStop):
            raise
        except Exception as exc:
            state['restore_errors'].append(f'{label}: {_mcm_describe(exc)}')

    if state['mcm_x_active']:
        yield from _try('mcm.x', bps.mv(mcm.x, cfg.mcm_xpos))
    if (state['hrm_state'] == HRM_IN and state['vertical_started']
            and state['mcm_ypos'] is not None):
        yield from _try('mcm.y', bps.mv(mcm.y, state['mcm_ypos']))
    elif state['mcm_y_return'] is not None:
        yield from _try('mcm.y', bps.mv(mcm.y, state['mcm_y_return']))
    if state['prep_started']:
        yield from mcm_setup_post(state)


def _mcm_result(state, status, reason):
    positions = {}
    for label, getter in (
            ('mcm.x', lambda: mcm.x.position), ('mcm.y', lambda: mcm.y.position),
            ('anc_xtal.y', lambda: anc_xtal.y.position),
            ('sample_stage.tx', lambda: sample_stage.tx.position),
            ('hrm2.dch', lambda: hrm2.dch.position),
            ('hrm2.dy', lambda: hrm2.dy.position)):
        try:
            positions[label] = getter()
        except Exception as exc:
            positions[label] = f'unreadable ({exc!r})'
    result = {
        'status': status, 'reason': reason, 'hrm_state': state['hrm_state'],
        'hrm_readbacks': state['hrm_readbacks'], 'axes': dict(state['axes']),
        'history': state['history'], 'corrections': state['corrections'],
        'scan_uids': state['uids'], 'final_positions': positions,
        'preparation_restored': state['restored'],
        'restore_errors': state['restore_errors'],
    }
    print('\n' + '*' * 60)
    print(f'MCM alignment: {status.upper()}' + (f' - {reason}' if reason else ''))
    for key in ('hrm_state', 'hrm_readbacks', 'axes', 'corrections',
                'final_positions', 'preparation_restored', 'restore_errors'):
        print(f'  {key}: {result[key]}')
    for rec in state['history']:
        print(f"  {rec['axis']} scan {rec['scan']}: max {rec['max']:.5f}, "
              f"error {rec['error']:+.5f}, uid {rec['uid']}")
    print('*' * 60 + '\n')
    return result


def mcm_setup(cfg=None, ask=input, **overrides):
    """Align the MCM horizontally and vertically (docs/mcm_alignment_scope.md).

    Usage: RE(mcm_setup()) or RE(mcm_setup(interactive=False, restore_positions=True)).
    Keyword overrides set McmAlignConfig fields. Returns a result dict.

    Failure (including a failed preparation move) restores the saved analyzer
    height, slits, anpd, whl and detector range. Operator interruption, abort
    or a declined prompt leaves everything as it is.
    """
    cfg = cfg or McmAlignConfig()
    for key, value in overrides.items():
        if not hasattr(cfg, key):
            raise TypeError(f'mcm_setup: unknown option {key!r}')
        setattr(cfg, key, value)
    state = _mcm_new_state(cfg)
    status, reason = 'success', ''
    try:
        yield from mcm_setup_prep(cfg, state)
        if _mcm_axis_enabled('horizontal', cfg.align_horizontal, cfg, ask):
            state['mcm_x_active'] = True
            yield from _mcm_axis_loop('horizontal', mcm.x, cfg.mcm_xpos,
                                      cfg.horizontal_tolerance, cfg, state, ask)
            state['axes']['horizontal'] = 'aligned'
        else:
            state['axes']['horizontal'] = 'skipped'
        if _mcm_axis_enabled('vertical', cfg.align_vertical, cfg, ask):
            state['vertical_started'] = True
            if state['hrm_state'] == HRM_OUT:
                yield from _mcm_vertical_hrm_out(cfg, state, ask)
            else:
                yield from _mcm_axis_loop('vertical', mcm.y, state['mcm_ypos'],
                                          cfg.vertical_tolerance, cfg, state, ask)
            state['axes']['vertical'] = 'aligned'
        else:
            state['axes']['vertical'] = 'skipped'
        if cfg.restore_positions:
            yield from mcm_setup_post(state)
    except (RequestAbort, RequestStop):
        raise  # operator interruption: no restoration
    except OperatorAbort as exc:
        status, reason = 'aborted', str(exc)
    except Exception as exc:
        status, reason = 'failed', _mcm_describe(exc)
        yield from _mcm_failure_cleanup(state)
    return _mcm_result(state, status, reason)


#*******************************************************************************************************
def analyzer_slit_scan(mtr, start, stop, gaps):
#    plt.clf()
    yield from bps.mv(mtr, 0)
    yield from bp.rel_scan([det2], mtr, start, stop, gaps)
    peak_stats = bec.peaks
    x0 = peak_stats['cen']['det2_current1_mean_value']
    print('*********************\n')
    print(f'{mtr.name} position center {x0}\n')
    if x0 > 1 or x0 < -1:
        print('*********************************************************\n')
        print(f'Verify the {mtr.name} data. Execution aborted!')
    else:
        yield from bps.mv(mtr, x0)
        mtr.set_current_position(0)


#*******************************************************************************************************
def san_setup():
#    acyy = ura_setup_prep()
    yield from analyzer_slit_scan(analyzer_slits.outboard, -1.2, 1.2, 41)
    yield from bps.mv(analyzer_slits.outboard, 2)

    yield from analyzer_slit_scan(analyzer_slits.inboard, -1.2, 1.2, 41)
    yield from bps.mv(analyzer_slits.inboard, -2)

    yield from analyzer_slit_scan(analyzer_slits.top, -1., 1., 41)
    yield from bps.mv(analyzer_slits.top, 2)
    
    yield from analyzer_slit_scan(analyzer_slits.bottom, -1., 1., 41)
    yield from bps.mv(analyzer_slits.bottom, -2)
    
#    ura_setup_post(acyy)
    print('*****************************************\n')
    print("Analyzer slits setup finished successfully\n")


#*******************************************************************************************************
def calculate_max_value(uid=-1, x="hrmE", y="lambda_det_stats7_total", delta=1, sampling=200):
    """
    This method gets a table (DataFrame) by using its uid. it finds the maximum value of the curve 
    under the sampled data by using the maximum y value and its neighboring data samples and then, 
    applying a polynomial regression over this curve. The model is used as an interpolation approach
    to generate more points between the original range and to return the x and y values of the
    maximum point of this new model

    Parameters
    ----------
    uid : int, optional
        id of the scan. The default is -1.
    x : str, optional
        label of the x values in the table. The default is "hrmE".
    channel : str, optional
        value of the channel with the y values. The default is 7.
    delta : int, optional
        total of points to be used on each side of the maximum value to generate the new model. The default is 1.
    sampling : int, optional
        total of sampling points to be used for interpolation. The default is 200.

    Raises
    ------
    ValueError
        The selected delta value is too big to be used based on the position of the maximum value in the table.

    Returns
    -------
    flaot
        x value of the maximum value.
    float
        y value of the maximum value.

    """

    run = trc[uid]
    table = run['primary'].read([x, y])
    #y = f'lambda_det_stats{channel}_total'
    
    #cp_df = df.copy()
    
    max_id = table[y].idxmax()
    
    # low limit check
    if max_id >= delta:
        low_max_id = max_id - delta
    else:
        raise ValueError("Delta value is greater than the lower limit of the dataset")
    
    # high limit check
    if max_id < len(table[y])-delta-1:
        high_max_id = max_id + delta + 1
    else:
        raise ValueError("Delta value is greater than the upper limit of the dataset")
    
    y_values = table[y][low_max_id:high_max_id]
    x_values = table[x][low_max_id:high_max_id]
    
    model = np.poly1d(np.polyfit(x_values, y_values, 2))
    
    resampled_x_values = np.linspace(x_values.iloc[0],x_values.iloc[-1],sampling)
    resampled_y_values = model(resampled_x_values)
    
    resample_df = pd.DataFrame({x:resampled_x_values, y:resampled_y_values})
    
    new_max_id = resample_df[y].idxmax()
    
    return resample_df[x][new_max_id], resample_df[y][new_max_id]


#*******************************************************************************************************
def LocalBumpSetup(silent=False):
    """
    Adjusts the e-beam local bump, i.e. horizontal & vertical positions of the x-ray beam on the XBPM1 screen
    Beam position at the XBPM2 on screen X=1245, Y=1066 with the CRL in

    """
    
    cond1 = strg_ring_orb_fb.uofb_pv.read()['srofb_uofb_pv']['value']
    cond2 = strg_ring_orb_fb.id_bump_pv.read()['srofb_id_bump_pv']['value']
    cond3 = strg_ring_orb_fb.nudge_pv.read()['srofb_nudge_pv']['value']
    pos3 = 0
    pos5 = 49
    cenX_target = 387
    cenY_target = 904
    max_dThe = 0.01  # mrad (= 10 urad)

    if cond1 != 2:
        print("****************** WARNING ******************")
        print("The UOFB is disabled. Operation is terminated")
        return
    if cond2 != 1:
        print("****************** WARNING ******************")
        print("The ID Bump is disabled. Operation is terminated")
        return
    if cond3 != 1:
        print("****************** WARNING ******************")
        print("The Nudge is disabled. Operation is terminated")
        return

    yield from bps.mv(bpm1_diag, pos5)
    cam1.cam.acquire_time.set(0.0002)
    cam1.cam.acquire_period.set(0.5)
    cam1.cam.num_images.set(1)
    cam1.cam.image_mode.set(2)
    cam1.cam.acquire.set(1)
    cam1.stats1.enable.set(1)
    cam1.stats1.compute_statistics.set(1)
    cam1.stats1.compute_centroid.set(1)
    Imax = cam1.stats1.max_value.get()
    if Imax < 2000:
        print('********************************')
        print('Low image intensity')
        print('Execution is terminated')
        return
    
    cam1.stats1.centroid_threshold.set(1000)
    centr = cam1.stats1.centroid.get()
    print(f'centroid X = {centr[1]:.2f}, target X = {cenX_target}\n')
    cenX = centr[1]
    dXc = cenX - cenX_target
    dThe = 1.e-3*dXc/6.0
    if abs(dThe) >= max_dThe:
        print("****************** WARNING ******************")
        print(f"Calculated horizontal e-beam shift {1.e3*dThe:0.1f} urad exceeds the limit of {1.e3*max_dThe:0.0f} urad. Horizontal correction skipped\n")
    else:
        print(f"Calculated horizontal e-beam shift {1.e3*dThe:0.1f} urad")
        if not silent:
            input_opts = input('Do you want to put it in (yes/no): ')
        else:
            input_opts = 'yes'
        
        if input_opts == 'yes':
            strg_ring_orb_fb.nudge_increment.set(dThe)
            strg_ring_orb_fb.horz_plane_nudge.set(1)
            yield from sleep(2)
            print('*****************************************')
            xtarget = strg_ring_orb_fb.xa_rbv.read()
            valx = 1.e3*xtarget['srofb_xa_rbv']['value']
            print(f'Horizontal angle correction was applied. New value {valx:0.1f} urad\n')
            print(strg_ring_orb_fb.nudge_status.alarm_status)
        else:
            print('*****************************************')
            print('Correction was canceled\n')

    crl_y_pos = crl.read()['crl_y']['value']
    if crl_y_pos < 1:
        print('\n')
        print('Error: CRL is in the x-ray beam. Vertical beam correction is canceled.\n')
        return
    
    print(f'centroid Y = {centr[0]:.2f}, target Y = {cenY_target}\n')
    cenY = centr[0]
    dYc = cenY - cenY_target
    dThe = -1.e-3*dYc/5.0
    if abs(dThe) >= max_dThe:
        print("****************** WARNING ******************")
        print(f"Calculated vertical e-beam shift {1.e3*dThe:0.1f} urad exceeds the limit of {1.e3*max_dThe:0.0f} urad. Vertical correction skipped\n")
    else:
        print(f"Calculated vertical e-beam shift {1.e3*dThe:0.1f} urad")
        if not silent:
            input_opts = input('Do you want to put it in (yes/no): ')
        else:
            input_opts = 'yes'
        
        if input_opts == 'yes':
            strg_ring_orb_fb.nudge_increment.set(dThe)
            strg_ring_orb_fb.vert_plane_nudge.set(1)
            yield from sleep(2)
            print('*****************************************')
            ytarget = strg_ring_orb_fb.xy_rbv.read()
            valy = 1.e3*ytarget['srofb_xy_rbv']['value']
            print(f'Vertical angle correction was applied. New value {valy:0.1f} urad\n')
            print(strg_ring_orb_fb.nudge_status.alarm_status)
        else:
            print('*****************************************')
            print('Correction was canceled\n')

    if not silent:
        update_opts = input('Do you want to move the XBPM1 back (yes/no): ')
    else:
        update_opts = 'yes'
    
    if update_opts == 'yes':
        cam1.cam.acquire.set(0)
        yield from bps.mv(bpm1_diag, pos3)


#*******************************************************************************************************
class CCRSetupError(RuntimeError):
    """Raised when the C crystal alignment cannot proceed or fails."""


def _ccr_fail(msg):
    print('*************************************')
    print(f'Error: {msg}. Execution aborted\n')
    raise CCRSetupError(msg)


def _ccr_check_hrm():
    # Must run before any motion: refuse to start if the HRM is in the beam
    rd = hrm2.read()
    if rd['hrm2_ux']['value'] > -5 or rd['hrm2_dx']['value'] > -5:
        _ccr_fail('HRM is in the beam')


def ccr_setup_prep():
# Prepares the URA for the C crystal setup
    yield from bps.mv(airpad, 1)
    yield from bps.mv(det2.em_range, 0)
    yield from bps.mv(spec.tth, 0)
    acyy = anc_xtal.read()['anc_xtal_y']['value']
    if acyy < 5:
        _ccr_fail('URA Y-position (acyy) is too low')

    yield from bps.mv(analyzer_slits.top, 0.1, analyzer_slits.bottom, -0.1, analyzer_slits.outboard, 1, analyzer_slits.inboard, -1)
    d21cnt = det2.current1.mean_value.read()['det2_current1_mean_value']['value']
    if d21cnt < 1.0e5:
        _ccr_fail('low intensity on D21')


#*******************************************************************************************************
def ccr_setup_post():
#   Recover positions after the C crystal alignment is finished
    yield from bps.mv(anpd, -90, analyzer_slits.top, 1, analyzer_slits.bottom, -1, analyzer_slits.outboard, 1.5, analyzer_slits.inboard, -1.5)


#*******************************************************************************************************
def ccr_the_setup():
#   Performs C crystal theta alignment; raises CCRSetupError if not converged
    x0 = 10
    kxmov = 0
    while abs(x0) > 5 and kxmov < 5:
        yield from bp.rel_scan([det2], analyzer.cfth, -150, 150, 31)
        x0 = bec.peaks['cen'].get('det2_current1_mean_value')
        if x0 is None:
            _ccr_fail('C crystal The-peak centre not found')
        # cfth is in urad; convention (sign included) confirmed on hardware
        x_deg = -np.rad2deg(1.e-6*x0)
        yield from bps.mvr(anc_xtal.the, x_deg)
        kxmov += 1

    if abs(x0) > 5:
        _ccr_fail('C crystal The-positioning did not converge')
    return kxmov


#*******************************************************************************************************
def _ccr_setup_body(s1, s2, s3):
    yield from ccr_setup_prep()

    # C crystal Theta alignment
    if s1 != 0:
        print('\n')
        print('Setting up the C crystal The-position\n')
        yield from bps.mv(det2.em_range, 0)
        yield from bps.sleep(1)
        yield from bps.mv(whl, 0, anpd, 40)
        yield from ccr_the_setup()

    # C crystal Y alignment
    if s2 != 0:
        print('\n')
        print('Setting up the C crystal Y-position\n')
        yield from bps.mv(det2.em_range, 1)
        yield from bps.sleep(1)
        yield from bps.mv(whl, 2, anpd, -145)
        yield from bps.mv(analyzer_slits.top, 0.04, analyzer_slits.bottom, -0.04)
        yield from bp.rel_scan([det2], anc_xtal.y, -0.3, 0.3, 61)
        peak_stats = bec.peaks
        p_max = peak_stats['max']['det2_current1_mean_value'][1]
        p_min = peak_stats['min']['det2_current1_mean_value'][1]
        if p_max - p_min < 1.e3:
            _ccr_fail('acyy position maximum not found')

        p_com = peak_stats['com']['det2_current1_mean_value']
        yield from bps.mv(anc_xtal.y, p_com)
        yield from bps.mv(analyzer_slits.top, 0.1, analyzer_slits.bottom, -0.1)

    # C crystal Chi alignment
    if s3 != 0:
        print('\n')
        print('Setting up the C crystal Chi-position\n')
        yield from bps.mv(whl, 0, anpd, -145)
        yield from bps.mv(det2.em_range, 0)
        yield from bps.sleep(1)
        yield from bps.mvr(anc_xtal.y, 0.6)
        yield from bps.mv(analyzer_slits.outboard, 0.5, analyzer_slits.inboard, -0.5)
        d22cnt = det2.current2.mean_value.read()['det2_current2_mean_value']['value']
        if d22cnt < 1.e5:
            _ccr_fail('low intensity on d22. Verify acchi manually')

        yield from bp.rel_scan([det2], analyzer.cchi, -0.15, 0.15, 31)
        yield from bps.mvr(analyzer.cchi, -0.15)
        x0 = bec.peaks['cen'].get('det2_current2_mean_value')
        if x0 is None:
            _ccr_fail('C crystal Chi-peak centre not found')
        yield from bps.mv(analyzer.cchi, x0)
        yield from bps.mv(analyzer_slits.outboard, 1, analyzer_slits.inboard, -1)

        # Re-checking C crystal The-position
        yield from bps.mv(anpd, 40)
        yield from ccr_the_setup()


#*******************************************************************************************************
def ccr_setup(s1=0, s2=0, s3=0):
#   C crystal alignment, namely the-position, analyzer vertical position and chi-position

    if s1 == 0 and s2 == 0 and s3 == 0:
        print("\n")
        print("Usage: ccr_setup(s1, s2, s3)")
        print("if s1 > 0, then execute acthe alignment, else - skip")
        print("if s2 > 0, then execute acyy alignment, else - skip")
        print("if s3 > 0, then execute cchi alignment, else - skip")
        return

    _ccr_check_hrm()  # before any motion, outside the restore wrapper
    yield from bpp.finalize_wrapper(_ccr_setup_body(s1, s2, s3), ccr_setup_post())
    print('\n')
    print('C crystal alignment finished\n')


#*******************************************************************************************************
def wcr_setup():
#   Performs W crystal alignment
    yield from bps.mv(anpd, -90, whl, 7, analyzer_slits.top, 0.1, analyzer_slits.bottom, -0.1, analyzer_slits.outboard, 1, analyzer_slits.inboard, -1)
    # yield from set_lambda_exposure(1)
    # yield from bp.rel_scan([lambda_det], analyzer.wfth, -20, 20, 41)
    yield from dscan(analyzer.wfth, -20, 20, 41, lambda_det, 1)
    x_pos = calculate_max_value(x="analyzer.wfth", sampling=100)
    yield from bps.mvr(analyzer.wfth, -20)
    yield from bps.mv(analyzer.wfth, x_pos)
    print('\n')
    print('W crystal alignment finished\n')


#*******************************************************************************************************
def hrm_in():
#   Moves HRM into the x-ray beam
    yield from bps.mv(hrm2.ux, 0, hrm2.dx, 0, hrm2.bs, 3)


#*******************************************************************************************************
def hrm_out():
#   Moves HRM out of the x-ray beam
    yield from bps.mv(hrm2.ux, -20, hrm2.dx, -20, hrm2.bs, 0)


#*******************************************************************************************************
def hrm_setup():
#   Performs HRM crystals alignment
    det4.em_range.set(0)
    det4.acquire_mode.set(0)
    det4.averaging_time.set(0.5)
    det5.em_range.set(0)
    det5.acquire_mode.set(0)
    det5.averaging_time.set(0.5)
    sleep(1)

    yield from hrm_in()
    yield from bps.mv(hrmE, 0, hrm2.d1, 0)
    yield from bps.mv(s1.top, 0.5, s1.bottom, -0.5, s1.outboard, 1, s1.inboard, -1)
    yield from bps.mv(hrm2.d2, 1, hrm2.d4, 0.7)
    yield from bps.mv(hrm2.d3, 2, hrm2.d5, 2)
    
    # 1st crystal alignment
    yield from bp.rel_scan([det4], hrm2.uth, -0.05, 0.05, 51)
    peak_stats = bec.peaks
    x_com = peak_stats['com']['det4_current2_mean_value']
    x_cen = peak_stats['cen']['det4_current2_mean_value']
    x_fwhm = peak_stats['fwhm']['det4_current2_mean_value']
    if abs(x_com-x_cen)/x_fwhm > 0.3:
        print('\n')
        print('Error: 1st crystal peak not found. Execution terminated.\n')
        return
    yield from bps.mv(hrm2.uth, x_cen, hrm2.d2, 0)

    # 2nd crystal alignment
    yield from bp.rel_scan([det4], hrm2.uif, -70, 70, 51)
    peak_stats = bec.peaks
    x_com = peak_stats['com']['det4_current3_mean_value']
    x_cen = peak_stats['cen']['det4_current3_mean_value']
    x_fwhm = peak_stats['fwhm']['det4_current3_mean_value']
    if abs(x_com-x_cen)/x_fwhm > 0.3:
        print('\n')
        print('Error: 2nd crystal peak not found. Execution terminated.\n')
        return
    yield from bps.mv(hrm2.uif, x_cen, hrm2.d3, 0)

    # 3rd crystal alignment
    yield from bp.rel_scan([det4], hrm2.dth, -0.01, 0.01, 51)
    peak_stats = bec.peaks
    x_com = peak_stats['com']['det4_current4_mean_value']
    x_cen = peak_stats['cen']['det4_current4_mean_value']
    x_fwhm = peak_stats['fwhm']['det4_current4_mean_value']
    if abs(x_com-x_cen)/x_fwhm > 0.3:
        print('\n')
        print('Error: 3rd crystal peak not found. Execution terminated.\n')
        return
    yield from bps.mv(hrm2.dth, x_cen, hrm2.d4, 0)

    # 4th crystal alignment
    yield from bp.rel_scan([det5], hrm2.dif, -70, 70, 51)
    peak_stats = bec.peaks
    x_com = peak_stats['com']['det5_current1_mean_value']
    x_cen = peak_stats['cen']['det5_current1_mean_value']
    x_fwhm = peak_stats['fwhm']['det5_current1_mean_value']
    if abs(x_com-x_cen)/x_fwhm > 0.3:
        print('\n')
        print('Error: 4th crystal peak not found. Execution terminated.\n')
        return
    yield from bps.mv(hrm2.dif, x_cen, hrm2.d5, 0)


#*******************************************************************************************************
def DxtalMesh(cnum=4, whl_pos=6, ctime=1):
    """
    Performs mesh scan of the analyzer D crystals: analyzer vertical position versus energy.
    
    Parameters
    -----------
    cnum    : int, optional
              number of measurements cycles
    whl_pos : int, optional
              position of the whl wheel
    ctime   : int, optional
              counting time for lambda detector
    
    """

    # spec_factory.prefix = "mesh_scan"
    yield from bps.mv(whl, whl_pos, spec.tth, 0)
    yield from bps.mv(analyzer_slits.top, 0.01, analyzer_slits.bottom, -0.01, analyzer_slits.outboard, 1.5, analyzer_slits.inboard, -1.5)
    yield from bps.mv(mcm_slits.outboard, 1.0, mcm_slits.inboard, -1.0)
    yield from set_lambda_exposure(ctime)
    acyy = anc_xtal.y.read()['anc_xtal_y']['value']
    hrmE_val = hrmE.read()['hrmE']['value']
#    print(acyy, hrmE_val)
    bec.enable_plots()
#    subs = [LiveGrid((150, 100), 'lambda_det_md7', xlabel='Energy', ylabel='acyy', ax=myfig.gca())]
#    plan = bpp.subs_wrapper(bp.rel_grid_scan([lambda_det], anc_xtal.y, -0.875, 0.625, 150, hrmE, -10, 10, 100, False), subs)
    for n in range(cnum):
        yield from bp.rel_grid_scan([lambda_det], anc_xtal.y, -0.875, 0.625, 150, hrmE, -10, 10, 100, False)
#        yield from plan
        sleep(600)
        
#    LiveGrid((150, 100), 'lambda_det_md7')
    bec.disable_plots()

#*******************************************************************************************************

def Beamline_Setup_1():
    # Performs 10-ID optics alignment down to the sample position
    #
    # 1. Move CRL to bypass position and wait for a few minutes for DCM to warm up.
    #
    print("Moving CRL out of the beam")
    cont_opts = input('Do you want to proceed (yes): ')
    if cont_opts == "" or cont_opts == "yes":
        yield from bps.mv(crl.y, 35.17)
        yield from sleep(2)
        print("CRL is out of the beam")
        cont_opts = input('Do you want to continue (yes): ')
        print(cont_opts)
        if cont_opts != "" and cont_opts != "yes":
            return
    #
    # 2. Check intensity at TM1 detector
    #
    counts = tm1.sum_all.mean_value.get()
    if counts < 1.0e-6:
        print("Low intensity at TM1 detector. Execution terminated")
        return
    #
    # 3. DCM setup
    #
    print("DCM setup")
    cont_opts = input('Do you want to proceed (yes): ')
    if cont_opts == "" or cont_opts == "yes":
        yield from dcm_setup()
        yield from sleep(2)
        print("DCM setup is done")
        cont_opts = input('Do you want to continue (yes): ')
        if cont_opts != "" and cont_opts != "yes":
            return
    #
    # 4. UGAP setup
    #
    print("Undulator gap setup")
    cont_opts = input('Do you want to proceed (yes): ')
    if cont_opts == "" or cont_opts == "yes":
        yield from ugap_setup()
        yield from sleep(2)
        print("Undulator gap setup is done")
        cont_opts = input('Do you want to continue (yes): ')
        if cont_opts != "" and cont_opts != "yes":
            return
    #
    # 5. Local bump setup
    #
    print("Local bump setup")
    cont_opts = input('Do you want to proceed (yes): ')
    if cont_opts == "" or cont_opts == "yes":
        yield from LocalBumpSetup(silent=True)
        yield from sleep(2)
        print("Local bump setup is done")
        cont_opts = input('Do you want to continue (yes): ')
        if cont_opts != "" and cont_opts != "yes":
            return
    #
    # 6. Move CRL to the focusing position and wait for a few minutes for DCM to warm up.
    #
    print("Moving CRL into the beam")
    cont_opts = input('Do you want to proceed (yes): ')
    if cont_opts == "" or cont_opts == "yes":
        yield from bps.mv(crl.y, 0.427)
        yield from sleep(2)
        print("CRL is in the beam")
        cont_opts = input('Do you want to continue (yes): ')
        if cont_opts != "" and cont_opts != "yes":
            return
    #
    # 7. DCM setup
    #
    print("DCM setup")
    cont_opts = input('Do you want to proceed (yes): ')
    if cont_opts == "" or cont_opts == "yes":
        yield from dcm_setup()
        yield from sleep(2)
        print("DCM setup is done")
        cont_opts = input('Do you want to continue (yes): ')
        if cont_opts != "" and cont_opts != "yes":
            return
        #
    # 8. CRL setup scan. Set CRL y-position to max intensity at TM1
    #
    print("CRL setup")
    cont_opts = input('Do you want to proceed (yes): ')
    if cont_opts == "" or cont_opts == "yes":
        # res = yield from dscan(crl.y, -0.2, 0.2, 20, tm1, 1, det_ch=[4])
        # fmax = res[0].max
        # xmax = fmax[0]
        yield from crl_setup()
        print("CRL setup is done")
        cont_opts = input('Do you want to continue (yes): ')
        if cont_opts != "" and cont_opts != "yes":
            return
    #
    # 9. Move HRM out of the beam
    #
    print("Moving HRM out of the beam")
    cont_opts = input('Do you want to proceed (yes): ')
    if cont_opts == "" or cont_opts == "yes":
        yield from hrm_out()
    
    print("Beamline setup 1 is completed successfully")