import bluesky.plans as bp
import bluesky.plan_stubs as bps
import bluesky.preprocessors as bpp
import bluesky.callbacks.fitting


def align_with_fit(dets, mtr, start, stop, gaps, md=None):
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

def set_lambda_exposure(exposure):
    det = lambda_det
    yield from bps.mv(det.cam.acquire_time, exposure, det.cam.acquire_period, exposure)

def check_zero(dets=None, start=-20, stop=20, gaps=200):
    ssxop = 0
    print('scanning zero')
    yield from bps.mov(spec.tth, 0, spec.phi, 0, sample_stage.sx, -1)
    sample_pos = yield from bps.read(sample_stage)
    print(sample_pos)
    if dets is None:
        dets = [det1]
    yield from bps.mov(whl, 7)
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
        yield from bps.mov(hrmE, target)

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

def some_plan(a, b ):
    ...