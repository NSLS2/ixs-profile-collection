import os
import bluesky.plans as bp
import bluesky.plan_stubs as bps
import bluesky.preprocessors as bpp
from bluesky.utils import RequestAbort, RequestStop
from ophyd import Device, EpicsSignal, EpicsSignalRO, Component as Cpt
from ophyd.areadetector import EpicsSignalWithRBV
import numpy as np
import pandas as pd
import lmfit
from bluesky.callbacks import LiveFit
from bluesky.callbacks.mpl_plotting import LiveGrid
from tabulate import tabulate
from utils.beamline_setup import (
    Step, beamline_setup_plan, print_step_table, beamline_setup_usage)
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
# def calc_stepup_fit(x):
#     # Calculates fitting parameters for step up function for MCM slits scan
#     y = 'det2_current1_mean_value'
#     run = trc[-1]
#     table = run['primary'].read([x, y])
#     model = lmfit.Model(stepup)
#     lf = LiveFit(model, y, {'x': x}, {'A': table[y].max(), 'sigma': 0.25, 'x0': 0, 'b':0})
#     for name, doc in run.documents():
#         lf(name, doc)
#     print(lf.result.values)
#     stup = stepup(table[x], **lf.result.values)
#     plt.clf()
#     plt.plot(table[x], table[y], label=f"raw data", marker = 'o', linestyle = 'none')
#     plt.plot(table[x], stup.values, label=f"data fit")
#     plt.legend()
#     return lf.result.values['x0']


#*******************************************************************************************************
# def calc_stepdwn_fit(x):
#     # Calculates fitting parameters for step down function for MCM slits scan
#     y = 'det2_current1_mean_value'
#     run = trc[-1]
#     table = run['primary'].read([x, y])
#     model = lmfit.Model(stepdown)
#     lf = LiveFit(model, y, {'x': x}, {'A': table[y].max(), 'sigma': 0.25, 'x0': 0, 'b':0})
#     for name, doc in run.documents():
#         lf(name, doc)
#     print(lf.result.values)
#     stdw = stepdown(table[x], **lf.result.values)
#     plt.clf()
#     plt.plot(table[x], table[y], label=f"raw data", marker = 'o', linestyle = 'none')
#     plt.plot(table[x], stdw.values, label=f"data fit")
#     plt.legend()
#     return lf.result.values['x0']


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
def _ugap_validate(stats, initial_gap):
    """Return xmax; raise AlignmentFailure with diagnostics."""
    MAX_SHIFT = 10.0        # um
    MIN_INTENSITY = 3.e-5
    fmax = stats.max
    diag = f'initial={initial_gap:.2f} max={fmax}'
    if fmax is None:
        raise AlignmentFailure(f'ID gap scan found no intensity maximum ({diag}).')
    xmax, ymax = fmax[0], fmax[1]
    if not ymax > MIN_INTENSITY:
        raise AlignmentFailure(
            f'ID gap maximum intensity {ymax:.3e} not above {MIN_INTENSITY:.1e} ({diag}).')
    if abs(xmax - initial_gap) > MAX_SHIFT:
        raise AlignmentFailure(
            f'ID gap maximum {xmax - initial_gap:+.2f} um from initial exceeds '
            f'+/-{MAX_SHIFT} um ({diag}).')
    return xmax


def ugap_setup(det=None, det_ch=(4,)):
    """Move the ivu22 gap to the intensity maximum. Raises AlignmentFailure on failure.

    det/det_ch select the scan detector and channels (default: tm1, channel 4).

    Failure restores the gap to its initial value and aborts the calling plan.
    """
    initial_gap = ivu22.position
    status, reason, restore_errors = 'success', '', []
    xmax = None
    try:
        res = yield from dscan(ivu22, -20, 20, 21, tm1 if det is None else det, 1,
                               det_ch=list(det_ch))
        # print("scan finished")
        xmax = _ugap_validate(res[0], initial_gap)
        print("Validation passed", xmax)
        yield from bps.mv(ivu22, xmax)
        print("Moved the gap to max")
    except (RequestAbort, RequestStop):
        raise
    except Exception as exc:
        status, reason = 'failed', _mcm_describe(exc)
        try:
            yield from bps.mv(ivu22, initial_gap)
        except (RequestAbort, RequestStop):
            raise
        except Exception as rexc:
            restore_errors.append(f'ivu22: {_mcm_describe(rexc)}')
    try:
        final = ivu22.position
    except Exception as exc:
        final = f'unreadable ({exc!r})'
    print('\n' + '*' * 60)
    print(f'ID gap setup: {status.upper()}' + (f' - {reason}' if reason else ''))
    print(f'  initial: {initial_gap}  final: {final}')
    print(f'  max position: {xmax}')
    print(f'  restore_errors: {restore_errors}')
    print('*' * 60 + '\n')
    if status == 'failed':
        raise AlignmentFailure(reason)


#*******************************************************************************************************
def _crl_validate(stats, initial_pos):
    """Return (xmax, contrast); raise AlignmentFailure with diagnostics."""
    MAX_SHIFT = 0.1
    MIN_CONTRAST = 0.08
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


CRL_IN_POS, CRL_OUT_POS = 0.427, 35.17
CRL_TOL = 0.01          # mm
CRL_IN_LIMIT = 1.0      # crl.y below this means the CRL is in the beam


def crl_in():
    """Move crl.y into the beam only if the CRL is out (crl.y >= 1 mm).

    If crl.y < 1 mm the position is left unchanged.
    Raises AlignmentFailure if the final position is not within +/-0.01 mm of the target.
    """
    pos = crl.y.position
    if pos < CRL_IN_LIMIT:
        print(f'CRL already in the beam (crl.y = {pos:.3f} mm). Position not changed.')
        return
    yield from bps.mv(crl.y, CRL_IN_POS)
    final = crl.y.position
    if abs(final - CRL_IN_POS) > CRL_TOL:
        raise AlignmentFailure(
            f'CRL in: crl.y = {final:.3f} mm, target {CRL_IN_POS} +/-{CRL_TOL} mm.')
    print(f'CRL moved in: crl.y {pos:.3f} -> {final:.3f} mm')


def crl_out():
    """Move crl.y out of the beam unless it is already at 35.17 +/-0.01 mm.

    Raises AlignmentFailure if the final position is not within +/-0.01 mm of the target.
    """
    pos = crl.y.position
    if abs(pos - CRL_OUT_POS) <= CRL_TOL:
        print(f'CRL already out of the beam (crl.y = {pos:.3f} mm). Position not changed.')
        return
    yield from bps.mv(crl.y, CRL_OUT_POS)
    final = crl.y.position
    if abs(final - CRL_OUT_POS) > CRL_TOL:
        raise AlignmentFailure(
            f'CRL out: crl.y = {final:.3f} mm, target {CRL_OUT_POS} +/-{CRL_TOL} mm.')
    print(f'CRL moved out: crl.y {pos:.3f} -> {final:.3f} mm')


def crl_setup():
    """Move crl.y to the intensity maximum. Raises AlignmentFailure on failure.

    Requires the CRL to be IN the beam (crl.y < 1 mm).
    Failure restores crl.y to its initial position and aborts the calling plan.
    """
    initial_pos = crl.y.position
    if not initial_pos < CRL_IN_LIMIT:
        raise AlignmentFailure(
            f'CRL setup requires the CRL IN the beam (crl.y < {CRL_IN_LIMIT} mm); '
            f'crl.y = {initial_pos:.3f} mm. Move it in with crl_in() first.')
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
    FWHM_MIN, FWHM_MAX = 44.0, 47.0
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
SAN_MIN_INTENSITY = 2.7e6
SAN_MIN_CONTRAST = 0.7
SAN_MAX_CENTER = 1.0
# blade: (open position, step direction, scan start, scan stop, points)
SAN_BLADES = {
    'outboard': (2, 'up', -1.2, 1.2, 41),
    'inboard': (-2, 'down', -1.2, 1.2, 41),
    'top': (2, 'up', -1.0, 1.0, 41),
    'bottom': (-2, 'down', -1.0, 1.0, 41),
}


def _san_step_fit(x, y, direction):
    """Fit a step edge with lmfit; return (x0, fit result)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    ymin, ymax = float(y.min()), float(y.max())
    span = float(x.max() - x.min())
    mid = 0.5 * (ymin + ymax)
    above = y > mid
    if direction == 'up':
        # stepup(): signal rises with x
        idx = int(np.argmax(above)) if above.any() else len(x) // 2
        model = lmfit.Model(stepup)
    else:
        # stepdown(): signal falls with x
        idx = int(np.argmax(~above)) if (~above).any() else len(x) // 2
        model = lmfit.Model(stepdown)
    params = model.make_params(A=ymax - ymin, sigma=0.1 * span, x0=float(x[idx]), b=ymin)
    params['sigma'].set(min=1e-6)
    result = model.fit(y, params, x=x)
    return float(result.params['x0'].value), result


def _san_validate(name, direction, x, y):
    """Fit the blade scan and validate it; return (x0, contrast)."""
    y = np.asarray(y, dtype=float)
    ymax, ymin = float(y.max()), float(y.min())
    contrast = (ymax - ymin) / (ymax + ymin) if (ymax + ymin) != 0 else None
    try:
        x0, result = _san_step_fit(x, y, direction)
    except Exception as exc:
        raise AlignmentFailure(f'{name}: step fit failed: {_mcm_describe(exc)}') from exc
    diag = f'x0={x0:.4f} contrast={contrast} max={ymax:.4g} min={ymin:.4g}'
    if not result.success:
        raise AlignmentFailure(f'{name}: step fit did not converge ({diag}).')
    if contrast is None or contrast < SAN_MIN_CONTRAST:
        raise AlignmentFailure(
            f'{name}: contrast {contrast} below {SAN_MIN_CONTRAST} ({diag}).')
    if abs(x0) > SAN_MAX_CENTER:
        raise AlignmentFailure(
            f'{name}: edge position {x0:.4f} outside +/-{SAN_MAX_CENTER} ({diag}).')
    return x0, contrast


def analyzer_slit_scan(mtr, direction, start, stop, gaps):
    """Scan one analyzer slit blade, fit the step edge, and center on it.

    Returns a record dict. Raises AlignmentFailure on failure; in that case
    the blade offset is not changed.
    """
    yield from bps.mv(mtr, 0)
    uid_holder = {}

    def _capture(name, doc):
        uid_holder['uid'] = doc['uid']

    plan = dscan(mtr, start, stop, gaps, det2, 1, det_ch=[0])
    yield from bpp.subs_wrapper(plan, {'start': [_capture]})
    uid = uid_holder.get('uid')
    if uid is None:
        raise AlignmentFailure(f'No Start document captured for the {mtr.name} scan.')
    try:
        flush_tiled()
        table = trc[uid]['primary'].read([mtr.name, 'det2_current1_mean_value'])
        x = table[mtr.name].values
        y = table['det2_current1_mean_value'].values
    except Exception as exc:
        raise AlignmentFailure(
            f'Reading {mtr.name} scan {uid} failed: {_mcm_describe(exc)}') from exc
    x0, contrast = _san_validate(mtr.name, direction, x, y)
    print('*********************\n')
    print(f'{mtr.name} edge position {x0:.4f}, contrast {contrast:.3f}\n')
    yield from bps.mv(mtr, x0)
    mtr.set_current_position(0)
    return {'x0': x0, 'contrast': contrast, 'uid': uid}


def _san_summary(status, reason, intensity, initial, records, restore_errors):
    final = {}
    for name in SAN_BLADES:
        try:
            final[name] = getattr(analyzer_slits, name).position
        except Exception as exc:
            final[name] = f'unreadable ({exc!r})'
    print('\n' + '*' * 60)
    print(f'Analyzer slits setup: {status.upper()}' + (f' - {reason}' if reason else ''))
    print(f'  open-slit intensity: {intensity}')
    print(f'  initial positions: {initial}')
    print(f'  final positions: {final}')
    for name, rec in records.items():
        print(f"  {name}: edge {rec['x0']:.4f}, contrast {rec['contrast']:.3f}, uid {rec['uid']}")
    print(f'  restore_errors: {restore_errors}')
    print('*' * 60 + '\n')


def san_setup():
    """Center the four analyzer slit blades on their step edges.

    The open-slit intensity is checked before any scan; if it fails the
    initial blade positions are restored. If a blade scan fails, blades
    already scanned stay where they are (offsets unchanged) and the failing
    blade is returned to its open position.
    """
    blades = {name: getattr(analyzer_slits, name) for name in SAN_BLADES}
    initial = {name: mtr.position for name, mtr in blades.items()}
    records, restore_errors = {}, []
    intensity = None
    status, reason = 'success', ''
    current = None
    try:
        # 1. Open all slits and check the intensity before any scan
        for name, (open_pos, *_rest) in SAN_BLADES.items():
            yield from bps.mv(blades[name], open_pos)
        intensity = det2.current1.mean_value.get()
        if not intensity > SAN_MIN_INTENSITY:
            current = 'initial'
            raise AlignmentFailure(
                f'Open-slit intensity {intensity:.4g} not above {SAN_MIN_INTENSITY:.1e}.')
        # 2. Scan the blades one by one
        for name, (open_pos, direction, start, stop, gaps) in SAN_BLADES.items():
            current = name
            records[name] = yield from analyzer_slit_scan(
                blades[name], direction, start, stop, gaps)
            yield from bps.mv(blades[name], open_pos)
            current = None
    except (RequestAbort, RequestStop):
        raise
    except Exception as exc:
        status, reason = 'failed', _mcm_describe(exc)
        if current == 'initial':
            # nothing scanned, nothing re-zeroed: restore initial positions
            for name, mtr in blades.items():
                try:
                    yield from bps.mv(mtr, initial[name])
                except (RequestAbort, RequestStop):
                    raise
                except Exception as rexc:
                    restore_errors.append(f'{name}: {_mcm_describe(rexc)}')
        elif current in SAN_BLADES:
            # failing blade only: back to its open position; scanned blades untouched
            try:
                yield from bps.mv(blades[current], SAN_BLADES[current][0])
            except (RequestAbort, RequestStop):
                raise
            except Exception as rexc:
                restore_errors.append(f'{current}: {_mcm_describe(rexc)}')
    _san_summary(status, reason, intensity, initial, records, restore_errors)
    if status == 'failed':
        raise AlignmentFailure(reason)
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
    # Tiled may return an xarray Dataset (no 'time' coordinate); the code
    # below needs a plain DataFrame with a RangeIndex.
    if hasattr(table, 'to_dataframe'):
        table = table.to_dataframe().reset_index(drop=True)
    else:
        table = pd.DataFrame(table).reset_index(drop=True)
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
# XBPM1 (BPM1 Y-foil) positions and e-beam bump settings
YAG_IN_POS, YAG_TOL = 49.0, 0.01        # mm, YAG screen in the beam
FOIL_OUT_POS, FOIL_OUT_TOL = 0.0, 0.1   # mm, foil out of the beam
BUMP_CENTROID_X, BUMP_CENTROID_Y = 387, 904
BUMP_CENTROID_TOL = 3.0                 # pixels, both axes
BUMP_MAX_ANGLE = 0.01                   # mrad (= 10 urad)
BUMP_MIN_IMAGE_MAX = 2000
CAM1_SETTLE_S = 1.0                     # let stats refresh after (re)starting
# (attribute path on cam1, value, tolerance)
CAM1_SETTINGS = (
    ('cam.acquire_time', 0.0002, 1e-6),
    ('cam.acquire_period', 0.5, 1e-3),
    ('cam.num_images', 1, 0),
    ('cam.image_mode', 2, 0),
    ('stats1.enable', 1, 0),
    ('stats1.compute_statistics', 1, 0),
    ('stats1.compute_centroid', 1, 0),
    ('stats1.centroid_threshold', 1000, 1e-6),
)


def _cam1_signal(path):
    obj = cam1
    for part in path.split('.'):
        obj = getattr(obj, part)
    return obj


def _cam1_apply_settings():
    for path, value, _ in CAM1_SETTINGS:
        yield from bps.mv(_cam1_signal(path), value)


def _cam1_verify_settings():
    """Read every camera setting back. Raise AlignmentFailure on any mismatch."""
    bad = []
    for path, value, tol in CAM1_SETTINGS:
        try:
            got = _cam1_signal(path).get()
            if abs(float(got) - value) > tol:
                bad.append(f'{path}={got} (expected {value})')
        except Exception as exc:
            bad.append(f'{path} unreadable ({exc!r})')
    if bad:
        raise AlignmentFailure('XBPM1 camera settings not as required: ' + '; '.join(bad))


def _cam1_stop():
    """Stop the XBPM1 camera acquisition and verify it by readback."""
    yield from bps.mv(cam1.cam.acquire, 0)
    if cam1.cam.acquire.get() != 0:
        raise AlignmentFailure('XBPM1 camera acquisition did not stop.')


def _xbpm1_yag_in():
    yield from bps.mv(bpm1_diag, YAG_IN_POS)
    pos = bpm1_diag.position
    if abs(pos - YAG_IN_POS) > YAG_TOL:
        raise AlignmentFailure(
            f'XBPM1 YAG position {pos:.4f} mm, required {YAG_IN_POS} +/-{YAG_TOL} mm.')


def _xbpm1_out():
    yield from bps.mv(bpm1_diag, FOIL_OUT_POS)
    pos = bpm1_diag.position
    if abs(pos - FOIL_OUT_POS) > FOIL_OUT_TOL:
        raise AlignmentFailure(
            f'XBPM1 foil position {pos:.4f} mm, required {FOIL_OUT_POS} +/-{FOIL_OUT_TOL} mm.')


def _bump_nudge(plane, d_angle, ask_text, silent, result):
    """Apply one e-beam angle correction (mrad). Raises if above the limit."""
    if abs(d_angle) >= BUMP_MAX_ANGLE:
        raise AlignmentFailure(
            f'Calculated {plane} e-beam shift {1.e3*d_angle:0.1f} urad exceeds the '
            f'limit of {1.e3*BUMP_MAX_ANGLE:0.0f} urad. Correction not applied.')
    print(f"Calculated {plane} e-beam shift {1.e3*d_angle:0.1f} urad")
    answer = 'yes' if silent else input(ask_text)
    if answer != 'yes':
        print('*****************************************')
        print('Correction was canceled\n')
        result['corrections'][plane] = 'canceled'
        return
    if plane == 'horizontal':
        pos = strg_ring_orb_fb.xa_rbv.read()['srofb_xa_rbv']['value']
        strg_ring_orb_fb.nudge_increment.set(d_angle)
        strg_ring_orb_fb.horz_plane_nudge.set(1)
    else:
        pos = strg_ring_orb_fb.xy_rbv.read()['srofb_xy_rbv']['value']
        strg_ring_orb_fb.nudge_increment.set(d_angle)
        strg_ring_orb_fb.vert_plane_nudge.set(1)
    yield from sleep(2)
    print('*****************************************')
    print(f'{plane.capitalize()} angle correction was applied. '
          f'Expected new value {1.e3*(pos + d_angle):0.1f} urad\n')
    print(strg_ring_orb_fb.nudge_status.alarm_status)
    result['corrections'][plane] = d_angle


def _local_bump_body(silent, keep_in_beam, state, result):
    yield from _xbpm1_yag_in()
    state['camera_started'] = True
    yield from _cam1_apply_settings()
    _cam1_verify_settings()          # before any ring correction
    cam1.cam.acquire.set(1)
    yield from bps.sleep(CAM1_SETTLE_S)

    Imax = cam1.stats1.max_value.get()
    if Imax < BUMP_MIN_IMAGE_MAX:
        raise AlignmentFailure(
            f'Low XBPM1 image intensity: max {Imax} < {BUMP_MIN_IMAGE_MAX}.')

    centr = cam1.stats1.centroid.get()
    cenX, cenY = centr[1], centr[0]
    result['centroid_before'] = {'x': cenX, 'y': cenY}
    print(f'centroid X = {cenX:.2f}, target X = {BUMP_CENTROID_X}')
    print(f'centroid Y = {cenY:.2f}, target Y = {BUMP_CENTROID_Y}\n')

    cond1 = strg_ring_orb_fb.uofb_pv.read()['srofb_uofb_pv']['value']
    cond2 = strg_ring_orb_fb.id_bump_pv.read()['srofb_id_bump_pv']['value']
    cond3 = strg_ring_orb_fb.nudge_pv.read()['srofb_nudge_pv']['value']
    problems = [text for ok, text in (
        (cond1 == 2, 'The UOFB is disabled'),
        (cond2 == 1, 'The ID Bump is disabled'),
        (cond3 == 1, 'The Nudge is disabled')) if not ok]
    if problems:
        within = (abs(cenX - BUMP_CENTROID_X) <= BUMP_CENTROID_TOL
                  and abs(cenY - BUMP_CENTROID_Y) <= BUMP_CENTROID_TOL)
        text = '; '.join(problems)
        if not within:
            raise AlignmentFailure(
                f'{text}, and the XBPM1 centroid ({cenX:.1f}, {cenY:.1f}) is not within '
                f'+/-{BUMP_CENTROID_TOL} px of the target '
                f'({BUMP_CENTROID_X}, {BUMP_CENTROID_Y}).')
        print('****************** WARNING ******************')
        print(f'{text}. Centroid is within +/-{BUMP_CENTROID_TOL} px of the target; '
              'no correction needed.\n')
        result['feedback_warning'] = text
        return

    yield from _bump_nudge('horizontal', 1.e-3*(cenX - BUMP_CENTROID_X)/6.0,
                           'Do you want to put it in (yes/no): ', silent, result)

    crl_y_pos = crl.read()['crl_y']['value']
    if crl_y_pos < 1:
        print('Note: CRL is in the x-ray beam. Vertical beam correction skipped.\n')
        result['corrections']['vertical'] = 'skipped (CRL in the beam)'
    else:
        yield from _bump_nudge('vertical', -1.e-3*(cenY - BUMP_CENTROID_Y)/5.0,
                               'Do you want to put it in (yes/no): ', silent, result)


def LocalBumpSetup(silent=False, keep_in_beam=False):
    """
    Adjusts the e-beam local bump, i.e. horizontal & vertical positions of the x-ray beam on the XBPM1 screen
    (centroid targets X=387, Y=904 px, +/-3 px).

    Raises AlignmentFailure for: YAG position error, camera setting mismatch, low
    image intensity, a correction above 10 urad, or disabled UOFB / ID bump / nudge
    when the centroid is not within +/-3 px of the target. A CRL in the beam skips
    only the vertical correction.

    silent       no prompts.
    keep_in_beam on success leave the YAG in and the camera acquiring (the caller
                 must take the image and stop the camera). On failure the camera is
                 always stopped and the foil moved out.
    Returns a dict with the centroid and the applied corrections.
    """
    state = {'camera_started': False, 'ok': False}
    result = {'corrections': {}}

    def cleanup():
        if state['ok'] and keep_in_beam:
            return
        if state['ok'] and not silent:
            if input('Do you want to move the XBPM1 back (yes/no): ') != 'yes':
                return
        for step_plan in ((_cam1_stop() if state['camera_started'] else None),
                          _xbpm1_out()):
            if step_plan is None:
                continue
            try:
                yield from step_plan
            except (RequestAbort, RequestStop):
                raise
            except Exception as exc:
                if state['ok']:
                    raise
                # do not mask the original failure
                print(f'LocalBumpSetup cleanup problem: {exc!r}')

    def body():
        yield from _local_bump_body(silent, keep_in_beam, state, result)
        state['ok'] = True

    yield from bpp.finalize_wrapper(body(), cleanup())
    return result


#*******************************************************************************************************
CCR_THE_FWHM = (110.0, 130.0)      # cfth units (urad)
CCR_THE_MIN_CONTRAST = 0.8
CCR_THE_MIN_INTENSITY = 3.2e6
CCR_Y_MIN_CONTRAST = 0.85
CCR_Y_MIN_INTENSITY = 3.0e6
CCR_CHI_MIN_CONTRAST = 0.8
CCR_CHI_MIN_INTENSITY = 1.0e6


def _ccr_fail(msg):
    print('*************************************')
    print(f'Error: {msg}. Execution aborted\n')
    raise AlignmentFailure(msg)


def _ccr_contrast(stats):
    try:
        ymax, ymin = stats.max[1], stats.min[1]
        if (ymax + ymin) != 0:
            return (ymax - ymin) / (ymax + ymin)
    except (TypeError, IndexError):
        pass
    return None


def _ccr_validate(label, stats, min_contrast, min_intensity, fwhm_range=None):
    """Return (cen, fwhm, contrast); raise AlignmentFailure with diagnostics."""
    cen, fwhm = stats.cen, getattr(stats, 'fwhm', None)
    contrast = _ccr_contrast(stats)
    ymax = None
    try:
        ymax = stats.max[1]
    except (TypeError, IndexError):
        pass
    diag = f'cen={cen} fwhm={fwhm} contrast={contrast} max={ymax}'
    if cen is None or ymax is None or contrast is None or (
            fwhm_range is not None and fwhm is None):
        _ccr_fail(f'{label} scan returned invalid results ({diag})')
    if not ymax > min_intensity:
        _ccr_fail(f'{label} maximum intensity {ymax:.3e} not above '
                  f'{min_intensity:.1e} ({diag})')
    if fwhm_range is not None and not fwhm_range[0] <= fwhm <= fwhm_range[1]:
        _ccr_fail(f'{label} FWHM {fwhm:.2f} outside '
                  f'[{fwhm_range[0]}, {fwhm_range[1]}] ({diag})')
    if not contrast > min_contrast:
        _ccr_fail(f'{label} contrast {contrast:.2f} not above '
                  f'{min_contrast} ({diag})')
    return cen, fwhm, contrast


def _ccr_the_validate(stats):
    return _ccr_validate('C crystal The', stats, CCR_THE_MIN_CONTRAST,
                         CCR_THE_MIN_INTENSITY, CCR_THE_FWHM)


def _ccr_y_validate(stats):
    return _ccr_validate('C crystal Y', stats, CCR_Y_MIN_CONTRAST,
                         CCR_Y_MIN_INTENSITY)


def _ccr_chi_validate(stats):
    return _ccr_validate('C crystal Chi', stats, CCR_CHI_MIN_CONTRAST,
                         CCR_CHI_MIN_INTENSITY)


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
#   Performs C crystal theta alignment; raises AlignmentFailure if not converged
    x0 = 10
    kxmov = 0
    while abs(x0) > 5 and kxmov < 5:
        res = yield from dscan(analyzer.cfth, -150, 150, 31, det2, 1, det_ch=[0])
        x0, _, _ = _ccr_the_validate(res[0])
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
        res = yield from dscan(anc_xtal.y, -0.3, 0.3, 61, det2, 1, det_ch=[0])
        p_cen, _, _ = _ccr_y_validate(res[0])
        yield from bps.mv(anc_xtal.y, p_cen)
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

        res = yield from dscan(analyzer.cchi, -0.15, 0.15, 31, det2, 1, det_ch=[1])
        yield from bps.mvr(analyzer.cchi, -0.15)
        x0, _, _ = _ccr_chi_validate(res[0])
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
WCR_FIELD = 'lambda_det_stats7_total'
WCR_MIN_INTENSITY = 40000
WCR_MIN_CONTRAST = 0.6
WCR_MAX_SHIFT = 5.0


def _wcr_validate(x, y, initial_pos, uid):
    """Fit the wfth scan and validate it; return (x_pos, y_fit, contrast, ymax)."""
    y = np.asarray(y, dtype=float)
    ymax, ymin = float(y.max()), float(y.min())
    contrast = (ymax - ymin) / (ymax + ymin) if (ymax + ymin) != 0 else None
    try:
        x_pos, y_fit = calculate_max_value(
            uid=uid, x='analyzer.wfth', y=WCR_FIELD, sampling=100)
        x_pos, y_fit = float(x_pos), float(y_fit)
    except Exception as exc:
        raise AlignmentFailure(
            f'Fit failed for analyzer.wfth scan {uid}: {_mcm_describe(exc)}') from exc
    diag = (f'initial={initial_pos:.3f} max_pos={x_pos:.3f} raw_max={ymax:.4g} '
            f'raw_min={ymin:.4g} contrast={contrast}')
    if not ymax > WCR_MIN_INTENSITY:
        raise AlignmentFailure(
            f'W crystal maximum intensity {ymax:.4g} not above {WCR_MIN_INTENSITY} ({diag}).')
    if contrast is None or contrast < WCR_MIN_CONTRAST:
        raise AlignmentFailure(
            f'W crystal contrast {contrast} below {WCR_MIN_CONTRAST} ({diag}).')
    if abs(x_pos - initial_pos) > WCR_MAX_SHIFT:
        raise AlignmentFailure(
            f'W crystal maximum {x_pos - initial_pos:+.3f} from initial exceeds '
            f'+/-{WCR_MAX_SHIFT} ({diag}).')
    return x_pos, y_fit, contrast, ymax


def wcr_setup():
    """Align the W crystal (analyzer.wfth) on the intensity maximum.

    Raises AlignmentFailure on failure and restores analyzer.wfth only.
    """
    initial_pos = analyzer.wfth.position
    status, reason, restore_errors = 'success', '', []
    x_pos = contrast = ymax = uid = None
    try:
        yield from bps.mv(anpd, -90, whl, 7, analyzer_slits.top, 0.1, analyzer_slits.bottom, -0.1, analyzer_slits.outboard, 1, analyzer_slits.inboard, -1)
        uid_holder = {}

        def _capture(name, doc):
            uid_holder['uid'] = doc['uid']

        plan = dscan(analyzer.wfth, -20, 20, 41, lambda_det, 1)
        yield from bpp.subs_wrapper(plan, {'start': [_capture]})
        uid = uid_holder.get('uid')
        if uid is None:
            raise AlignmentFailure('No Start document captured for the analyzer.wfth scan.')
        try:
            flush_tiled()
            table = trc[uid]['primary'].read(['analyzer.wfth', WCR_FIELD])
            x = table['analyzer.wfth'].values
            y = table[WCR_FIELD].values
        except Exception as exc:
            raise AlignmentFailure(
                f'Reading analyzer.wfth scan {uid} failed: {_mcm_describe(exc)}') from exc
        x_pos, _y_fit, contrast, ymax = _wcr_validate(x, y, initial_pos, uid)
        yield from bps.mvr(analyzer.wfth, -20)
        yield from bps.mv(analyzer.wfth, x_pos)
    except (RequestAbort, RequestStop):
        raise
    except Exception as exc:
        status, reason = 'failed', _mcm_describe(exc)
        try:
            yield from bps.mv(analyzer.wfth, initial_pos)
        except (RequestAbort, RequestStop):
            raise
        except Exception as rexc:
            restore_errors.append(f'analyzer.wfth: {_mcm_describe(rexc)}')
    try:
        final = analyzer.wfth.position
    except Exception as exc:
        final = f'unreadable ({exc!r})'
    print('\n' + '*' * 60)
    print(f'W crystal setup: {status.upper()}' + (f' - {reason}' if reason else ''))
    print(f'  initial: {initial_pos}  final: {final}')
    print(f'  max position: {x_pos}  raw max: {ymax}  contrast: {contrast}  uid: {uid}')
    print(f'  restore_errors: {restore_errors}')
    print('*' * 60 + '\n')
    if status == 'failed':
        raise AlignmentFailure(reason)
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
def _hrm_peak(label, stats):
    """Return the peak centre; raise AlignmentFailure if no valid peak was found."""
    com, cen, fwhm = stats.com, stats.cen, stats.fwhm
    diag = f'com={com} cen={cen} fwhm={fwhm}'
    if com is None or cen is None or fwhm is None or not fwhm > 0:
        raise AlignmentFailure(f'HRM {label} crystal: invalid peak statistics ({diag}).')
    if abs(com - cen)/fwhm > 0.3:
        raise AlignmentFailure(f'HRM {label} crystal peak not found ({diag}).')
    return cen


def _hrm_crystal_scan(label, motor, start, stop, det, channel, d_motor, d_value):
    """Scan one HRM crystal angle with dscan, move to the peak, then set the D motor."""
    res = yield from dscan(motor, start, stop, 51, det, 0.5, det_ch=[channel])
    cen = _hrm_peak(label, res[0])
    yield from bps.mv(motor, cen, d_motor, d_value)


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
    
    # 1st..4th crystal alignment. det_ch indexes det.hints['fields']
    # (det4: current1..4 -> 0..3; det5: current1 -> 0).
    yield from _hrm_crystal_scan('1st', hrm2.uth, -0.05, 0.05, det4, 1, hrm2.d2, 0)
    yield from _hrm_crystal_scan('2nd', hrm2.uif, -70, 70, det4, 2, hrm2.d3, 0)
    yield from _hrm_crystal_scan('3rd', hrm2.dth, -0.01, 0.01, det4, 3, hrm2.d4, 0)
    yield from _hrm_crystal_scan('4th', hrm2.dif, -70, 70, det5, 0, hrm2.d5, 0)


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

#*******************************************************************************************************
# Beamline_Setup: production sequence run by the framework in utils/beamline_setup.py.
# Not validated on hardware. Steps not yet implemented (lambda tdy, HRM/MCM matching,
# ID gap re-optimisation with det5, energy resolution check) are NOT in the table.
BEAMLINE_SETUP_REPORT_DIR = '/nsls2/data/ixs/legacy/beamline_setup_reports'
BPM1_TIFF_PREFIX = 'XF:10IDA-BI{BPM1-Cam:1}TIFF1:'
BPM1_IMAGE_DIR = '/nsls2/data/ixs/legacy/Xbpm1/'
BPM1_IMAGE_NAME = 'bpm1'
BPM1_IMAGE_TEMPLATE = '%s%s_%3.3d.tif'
BPM1_FRAME_WAIT_S = 1.0     # frames reach the TIFF plugin only after it is enabled
BPM1_FILE_WAIT_S = 5.0


class Bpm1Tiff(Device):
    """TIFF file plugin of the XBPM1 camera (created only when the image is taken)."""
    enable = Cpt(EpicsSignalWithRBV, 'EnableCallbacks')
    file_path = Cpt(EpicsSignalWithRBV, 'FilePath', string=True)
    file_name = Cpt(EpicsSignalWithRBV, 'FileName', string=True)
    file_template = Cpt(EpicsSignalWithRBV, 'FileTemplate', string=True)
    write_mode = Cpt(EpicsSignalWithRBV, 'FileWriteMode')
    auto_save = Cpt(EpicsSignalWithRBV, 'AutoSave')
    auto_increment = Cpt(EpicsSignalWithRBV, 'AutoIncrement')
    write_file = Cpt(EpicsSignal, 'WriteFile')
    full_name = Cpt(EpicsSignalRO, 'FullFileName_RBV', string=True)
    write_message = Cpt(EpicsSignalRO, 'WriteMessage', string=True)
    unique_id = Cpt(EpicsSignalRO, 'UniqueId_RBV')


def _bpm1_tiff_setup(tiff):
    """Apply and verify the TIFF settings. Raises RuntimeError on a mismatch."""
    yield from bps.mv(tiff.file_path, BPM1_IMAGE_DIR, tiff.file_name, BPM1_IMAGE_NAME,
                      tiff.file_template, BPM1_IMAGE_TEMPLATE, tiff.write_mode, 0,
                      tiff.auto_save, 0, tiff.auto_increment, 1)
    bad = []
    if str(tiff.file_path.get()).rstrip('/') != BPM1_IMAGE_DIR.rstrip('/'):
        bad.append(f'file_path={tiff.file_path.get()!r}')
    if tiff.file_name.get() != BPM1_IMAGE_NAME:
        bad.append(f'file_name={tiff.file_name.get()!r}')
    if tiff.file_template.get() != BPM1_IMAGE_TEMPLATE:
        bad.append(f'file_template={tiff.file_template.get()!r}')
    for label, sig, want in (('write_mode', tiff.write_mode, 0),
                             ('auto_save', tiff.auto_save, 0),
                             ('auto_increment', tiff.auto_increment, 1)):
        if int(sig.get()) != want:
            bad.append(f'{label}={sig.get()} (expected {want})')
    if bad:
        raise RuntimeError('TIFF settings not as required: ' + '; '.join(bad))


def bpm1_orbit_image(report_dir=BEAMLINE_SETUP_REPORT_DIR):
    """Read back the orbit and save the last XBPM1 camera frame (camera still acquiring).

    Raises AlignmentFailure if the camera settings are not as required. A failure
    to save the image only records image_path=None and an image_note. The TIFF
    callbacks are restored and the camera acquisition is always stopped afterwards.
    Returns a JSON-safe dict for the report.
    """
    state = {'tiff': None, 'enable_prev': None, 'camera_touched': True, 'ok': False}
    result = {}

    def body():
        result['orbit_x_angle_rbv'] = float(
            strg_ring_orb_fb.xa_rbv.read()['srofb_xa_rbv']['value'])
        result['orbit_y_angle_rbv'] = float(
            strg_ring_orb_fb.xy_rbv.read()['srofb_xy_rbv']['value'])
        result['uofb'] = int(strg_ring_orb_fb.uofb_pv.get())
        result['id_bump'] = int(strg_ring_orb_fb.id_bump_pv.get())
        result['nudge'] = int(strg_ring_orb_fb.nudge_pv.get())
        result['nudge_status'] = str(strg_ring_orb_fb.nudge_status.get())
        _cam1_verify_settings()
        centr = cam1.stats1.centroid.get()
        cx, cy = float(centr[1]), float(centr[0])
        result.update(
            centroid_x=cx, centroid_y=cy,
            target_x=BUMP_CENTROID_X, target_y=BUMP_CENTROID_Y,
            centroid_deviation={'x': cx - BUMP_CENTROID_X, 'y': cy - BUMP_CENTROID_Y},
            image_max=float(cam1.stats1.max_value.get()),
            image_path=None, image_path_relative=None, image_note='')
        try:
            tiff = Bpm1Tiff(BPM1_TIFF_PREFIX, name='bpm1_tiff')
            tiff.wait_for_connection(timeout=5)
            state['tiff'] = tiff
            yield from _bpm1_tiff_setup(tiff)
            state['enable_prev'] = tiff.enable.get()
            yield from bps.mv(tiff.enable, 1)
            yield from bps.sleep(BPM1_FRAME_WAIT_S)
            try:
                result['frame_unique_id'] = int(tiff.unique_id.get())
            except Exception:
                result['frame_unique_id'] = None
            yield from bps.mv(tiff.write_file, 1)
            full = str(tiff.full_name.get())
            if not full:
                raise RuntimeError('IOC reported no file name '
                                   f'(WriteMessage: {tiff.write_message.get()!r})')
            result['image_path'] = full
            result['image_path_relative'] = os.path.relpath(full, report_dir)
            for _ in range(int(BPM1_FILE_WAIT_S / 0.5) + 1):
                if os.path.exists(full):
                    break
                yield from bps.sleep(0.5)
            else:
                result['image_note'] = 'file not visible from this host'
        except (RequestAbort, RequestStop):
            raise
        except Exception as exc:
            result['image_path'] = result['image_path_relative'] = None
            result['image_note'] = f'image not saved: {exc!r}'
            print(f'WARNING: XBPM1 {result["image_note"]}')
        state['ok'] = True

    def cleanup():
        problems = []
        tiff = state['tiff']
        if tiff is not None and state['enable_prev'] is not None:
            try:
                yield from bps.mv(tiff.enable, state['enable_prev'])
            except (RequestAbort, RequestStop):
                raise
            except Exception as exc:
                problems.append(f'TIFF callbacks not restored: {exc!r}')
        try:
            yield from _cam1_stop()
            result['camera_stopped'] = True
        except (RequestAbort, RequestStop):
            raise
        except Exception as exc:
            result['camera_stopped'] = False
            problems.append(f'camera not stopped: {exc!r}')
        if problems:
            if state['ok']:
                raise AlignmentFailure('; '.join(problems))
            print('bpm1_orbit_image cleanup problem: ' + '; '.join(problems))

    yield from bpp.finalize_wrapper(body(), cleanup())
    return result


def _no_prompt(*args, **kwargs):
    raise RuntimeError('prompt requested in silent mode')


def _alignment_scan_step(keyword):
    """alignment_scan wrapper: its precondition failures become AlignmentFailure."""
    def factory():
        failed = []
        for check in _SCAN_REGISTRY[keyword]['checks']:
            ok, msg = check()
            if not ok:
                failed.append(msg)
        if failed:
            raise AlignmentFailure(f'alignment_scan({keyword!r}) conditions not met: '
                                   + ' '.join(failed))
        yield from alignment_scan(keyword)
    return factory


def beamline_setup_step_table(report_dir=BEAMLINE_SETUP_REPORT_DIR, interactive=True,
                              dcm_warmup_s=150, crl_warmup_s=150,
                              orbit_settle_s=120, hrm_warmup_s=1800):
    """Return the production step table (list of Step). Moves nothing."""
    run = {'crl_was_out': None}

    def crl_out_step():
        run['crl_was_out'] = abs(crl.y.position - CRL_OUT_POS) <= CRL_TOL
        yield from crl_out()

    def wait(seconds, label):
        def factory():
            print(f'{label}: waiting {seconds} s')
            yield from bps.sleep(seconds)
            return {'waited_s': seconds}
        return factory

    def dcm_warmup():
        seconds = 0 if run['crl_was_out'] is True else dcm_warmup_s
        print(f'DCM warm-up: waiting {seconds} s '
              f'(CRL already out: {run["crl_was_out"]})')
        yield from bps.sleep(seconds)
        return {'waited_s': seconds, 'crl_was_out': run['crl_was_out']}

    def bump():
        res = yield from LocalBumpSetup(silent=not interactive, keep_in_beam=True)
        return res

    def image():
        res = yield from bpm1_orbit_image(report_dir)
        return res

    def mcm_step():
        ask = input if interactive else _no_prompt
        res = yield from mcm_setup(interactive=interactive, ask=ask)
        if res['status'] == 'failed':
            raise AlignmentFailure(res['reason'])
        if res['status'] == 'aborted':
            raise OperatorAbort(res['reason'])
        return {'status': res['status'], 'axes': res['axes'],
                'final_positions': res['final_positions']}

    def ccr_step():
        yield from ccr_setup(1, 1, 0)

    entries = [
        ('XBPM1 foil out [bpm1_diag -> 0]', 'auto', _xbpm1_out),
        ('HRM out [hrm_out]', 'auto', hrm_out),
        ('CRL out [crl_out]', 'auto', crl_out_step),
        ('DCM warm-up wait', 'auto', dcm_warmup),
        ('DCM setup [dcm_setup]', 'auto', dcm_setup),
        ('ID gap setup [ugap_setup]', 'auto', ugap_setup),
        ('DCM alignment scan', 'optional', _alignment_scan_step('dcm')),
        ('E-beam local bump, XBPM1 on YAG [LocalBumpSetup]', 'auto', bump),
        ('Orbit settle wait', 'auto', wait(orbit_settle_s, 'Orbit settle')),
        ('XBPM1 orbit readback and image', 'auto', image),
        ('XBPM1 foil out [bpm1_diag -> 0]', 'auto', _xbpm1_out),
        ('CRL in [crl_in]', 'auto', crl_in),
        ('CRL warm-up wait', 'auto', wait(crl_warmup_s, 'CRL warm-up')),
        ('DCM setup [dcm_setup]', 'auto', dcm_setup),
        ('CRL y setup [crl_setup]', 'auto', crl_setup),
        ('DCM alignment scan (CRL in)', 'optional', _alignment_scan_step('dcm')),
        ('MCM alignment [mcm_setup]', 'auto', mcm_step),
        ('MCM alignment scan', 'optional', _alignment_scan_step('mcm')),
        ('Analyzer slits setup [san_setup]', 'auto', san_setup),
        ('C crystal setup [ccr_setup 1 1 0]', 'auto', ccr_step),
        ('CCR alignment scan', 'optional', _alignment_scan_step('ccr')),
        ('W crystal setup [wcr_setup]', 'auto', wcr_setup),
        ('WCR alignment scan', 'optional', _alignment_scan_step('wcr')),
        ('HRM in [hrm_in]', 'auto', hrm_in),
        ('HRM warm-up wait', 'auto', wait(hrm_warmup_s, 'HRM warm-up')),
        ('HRM setup [hrm_setup]', 'auto', hrm_setup),
        ('HRM alignment scan', 'optional', _alignment_scan_step('hrm')),
    ]
    table, n = [], 0
    for name, kind, factory in entries:
        if kind == 'optional':
            table.append(Step(f'{n}a', name, kind, factory))
        else:
            n += 1
            table.append(Step(str(n), name, kind, factory))
    return table


def Beamline_Setup(start=1, stop=None, interactive=True, optional_scans=False,
                   confirm=input, report_dir=BEAMLINE_SETUP_REPORT_DIR,
                   first_run=True, dcm_warmup_s=150, crl_warmup_s=150,
                   orbit_settle_s=120, hrm_warmup_s=1800):
    """Automatic beamline alignment. Usage: RE(Beamline_Setup()).

    first_run=True (default) requires interactive=True: every step is confirmed.
    start/stop only select a range of steps (not a validated resume).
    See Beamline_Setup_usage() and Beamline_Setup_table().
    """
    if first_run and not interactive:
        raise ValueError('first_run=True requires interactive=True '
                         '(pass first_run=False for a silent run).')
    waits = (dcm_warmup_s, crl_warmup_s, orbit_settle_s, hrm_warmup_s)
    if not all(isinstance(w, (int, float)) and not isinstance(w, bool) and w >= 0
               for w in waits):
        raise ValueError('wait durations must be numbers >= 0')
    table = beamline_setup_step_table(report_dir, interactive, dcm_warmup_s,
                                      crl_warmup_s, orbit_settle_s, hrm_warmup_s)
    return beamline_setup_plan(
        table, start=start, stop=stop, interactive=interactive,
        optional_scans=optional_scans, confirm=confirm, report_dir=report_dir,
        failure_types=(AlignmentFailure,), abort_types=(OperatorAbort,))


def Beamline_Setup_table(start=1, stop=None, optional_scans=False):
    """Print the Beamline_Setup step table. Runs nothing."""
    return print_step_table(beamline_setup_step_table(), start, stop, optional_scans)


def Beamline_Setup_usage():
    """Print how to use Beamline_Setup."""
    return beamline_setup_usage()
