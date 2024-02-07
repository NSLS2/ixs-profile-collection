import bluesky.plans as bp
import bluesky.plan_stubs as bps
import bluesky.preprocessors as bpp
import bluesky.callbacks.fitting
import numpy as np
import pandas as pd
import lmfit
from bluesky.callbacks import LiveFit
from bluesky.suspenders import SuspendFloor
from ophyd import EpicsSignal
from tabulate import tabulate

# tm1sum = EpicsSignal('XF:10ID-BI:TM176:SumAll:MeanValue_RBV')
# susp = SuspendFloor(tm1sum, 1.e-5, resume_thresh = 1.e-5, sleep = 1*60)

# uofb_pv = EpicsSignal("SR:UOFB{}ConfigMode-I", name="uofb_pv")
# id_bump_pv = EpicsSignal("SR:UOFB{C10-ID}Enabled-I", name="id_bump_pv")
# nudge_pv = EpicsSignal("SR:UOFB{C10-ID}Nudge-Enabled", name="nudge_pv")


def align_with_fit(dets, mtr, start, stop, gaps, md=None):
    # Performs relative scan of motor and retuns data staistics

    md = md or {}

    local_peaks = []
    for det in dets:
        for hint in det.hints['fields']:
            local_peaks.append(
                bluesky.callbacks.fitting.PeakStats(mtr.hints['fields'][0], hint)
            ) 
    # TODO use relative wrapper to avoid the reset behavior (or make it optional)
    plan = bpp.subs_wrapper(
        bp.rel_scan(dets, mtr, start, stop, gaps+1, md=md), 
        local_peaks
        )
    yield from plan
    return local_peaks

#def set_lambda_exposure(exposure):
#    det = lambda_det
#    yield from bps.mv(det.cam.acquire_time, exposure, det.cam.acquire_period, exposure)

def check_zero(dets=None, start=-20, stop=20, gaps=200, exp_time=1):
    # Performs relative scan of the HRM energy at tth = 0 and positions it to the peak center

    #ssxop = 0
    print('scanning zero')
    #yield from bps.mov(spec.tth, 0, spec.phi, 0, sample_stage.sx, -1)
    yield from bps.mv(spec.tth, 0)
    sample_pos = yield from bps.read(sample_stage)
    print(sample_pos)
    if dets is None:
        dets = [lambda_det]
        yield from set_lambda_exposure(exp_time)

    yield from bps.mv(whl, 7)
    for d in dets:
        # set the exposure times
        pass

    local_peaks = yield from align_with_fit(dets, hrmE, start, stop, gaps)
    
    cen = local_peaks[0].cen
    if cen is not None:
        target = 0.2 * (cen // .2)
        # move too far for backlash compensation
        yield from bps.mv(hrmE, target - 20)
        # apporach target from negative side 
        yield from bps.mv(hrmE, target)

def do_the_right_thing(i_time):
    yield from bps.mv(det1.integration_time, i_time)
    yield from count([det1])

def ct(exp_time):
    yield from bps.mv(sclr.preset_time, exp_time)
    yield from bp.count([sclr])


def double_ct(exp_time):
    yield from ct(exp_time)
    # yield from bps.mv(sample_stage.sx, 0)
    yield from ct(exp_time)

def Lipid_Qscan():
    # Test plan for the energy scan at several Q values
    tth001 = 16.8
    Qq = [1, 2, 3]
    c22 = sclr.channels.chan22
    yield from bps.mv(analyzer_slits.top, 1, analyzer_slits.bottom, -1, analyzer_slits.outboard, 1.5, analyzer_slits.inboard, -1.5)
    for kk in range(6):
        yield from bps.mv(anapd, 25)
        #yield from set_lambda_exposure(2)
        yield from check_zero(exp_time=2)
        yield from bps.mv(whl, 0)

        for q in Qq:
            th = qq2th(q)
            yield from bps.mv(spec.tth, th)
            yield from hrmE_dscan(-15, 15, 150, 30)

            yield from bps.mvr(sample_stage.sx, 0.03)
            yield from bps.mv(spec.tth, tth001)
            yield from set_lambda_exposure(5)
            loc_peaks = yield from align_with_fit([lambda_det], sample_stage.sy, -0.1, 0.1, 40)
            max_pos = local_peaks[0].max
            yield from bps.mvr(sample_stage.sy, -0.1)
            yield from bps.mv(sample_stage.sy, max_pos)
            loc_peaks = yield from align_with_fit([lambda_det], sample_stage.sz, -2, 2, 40)
            max_pos = local_peaks[0].max
            yield from bps.mv(sample_stage.sz, max_pos)
            loc_peaks = yield from align_with_fit([lambda_det], sample_stage.sy, -0.1, 0.1, 40)
            max_pos = local_peaks[0].max
            yield from bps.mvr(sample_stage.sy, -0.1)
            yield from bps.mv(sample_stage.sy, max_pos)

        yield from bps.mv(anapd, 3, spec.tth, 1)
        yield from bps.mv(sclr.channels.chan22.preset_time, 5)
        yield from bp.scan([c22], spec.tth, 1, 21, 101)

def Lipid_Qscan_wBC():
    # Lipid_Qscan with beam check
    yield from bpp.suspend_wrapper(Lipid_Qscan(), susp)



