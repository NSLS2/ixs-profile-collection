from bluesky.plans import scan
from bluesky.callbacks.fitting import PeakStats
from bluesky.preprocessors import subs_wrapper

def check_zero():
    ssxop = 0
    E0
    print('scanning zero')
    yield from bps.mov(tth, 0, phi, 0, ssx, -1)
    sample_pos = yield from bps.read(sample_stage)
    print(sample_pos)
    dets = [det1]
    yield from bps.mov(whl, 5)
    yield from bps.mov(hrmE, E0 - 20)
    for d in dets:
        # set the exposure times
        pass

    local_peaks = bluesky.callbacks.fitting.PeakStats(
        hrmE.energy.name, dets[0].hints['fields'][0])
    plan = bpp.subs_wrapper(rel_scan(dets, hrmE, 0, 40, 200,
                            md={'reason': 'alignment'}), [local_peaks])

    yield from plan
    cen = local_peaks.cen
    target = 0.2 * (cen // .2)
    yield from bps.mov(hrmE, target)


def generic_alignment_plan(det, motor, start, stop, num_steps,
                           move_centre=True, move_tolerance=0.2):
    '''A generic alignment plan that scans 'motor' and looks for a peak centre.

    This scans motor between 'start' and 'stop' with 'num_steps' and at the end
    looks for a peak using the ``bluesky.calbacks.fitting.PeakStats`` callback.
    It will then move the motor to the 'peak' and 'zero' the axis provided
    'move_centre' is True and the new 'centre' is within 'move_tolerance' of
    the starting location.

    It is expected that this will be used with ``functools.partial`` in order
    to generate a number of different 'alignment' plans using the following:

    .. code::

        from functools import partial
        align_motor_x = partial(
            generic_alignment_plan(some_detector, motor_x, some_start,
                                   some_stop, some_num_steps,
                                   move_tolerance=some_move_tolerance))

    This motor alignment scan can then be called using

    .. code::

        RE(align_motor_x())
        # or
        RE(align_motor_x(move_centre=False))

    Parameters
    ----------
    det : ophyd_obj
        The ophyd object to 'read' at each step in the scan
    motor : ophyd_obj
        The ophyd object to 'scan'
    start : float
        The start value for the scan
    stop : float
        The stop value for the scan
    num_steps : int
        The number of steps in the scan
    move_centre : bool, optional
        Indicates if the plan should move to the centre at the end, default is
        True.
    move_tolerance : float, optional
        The tolerance vlaue (i.e. if the centre is not within 'move_float' of
        the current location do not move but instead raise an error).
    '''
    # setup a local version of PeakStats to find the peak centre
    local_peaks = PeakStats(motor.name, det.hints['fields'][0])
    # define the plan.
    plan = subs_wrapper(scan([det], motor, start, stop, num_steps,
                             md={'reason': 'alignment'}),
                        [local_peaks])
    # find the initial location of motor
    initial_pos = motor.read()[motor.name]['value']

    # yield from each of the part of the plan
    uid = yield from plan
    if uid:  # This ensures that ``summarize_plan`` will work on this plan
        target = num_steps * (local_peaks.cen // num_steps)  # nearest step
        if (move_centre and abs(target-initial_pos)<=move_tolerance):
            yield from bps.mov(motor, target)
    elif move_centre:  # if move_centre add a print statemnt for summarize_plan
        print (f'may move {motor.name} to the centre value if it is within '
               f'{move_tolerance} of the initial motor value')
