from bluesky.utils import make_decorator
from bluesky.plan_stubs import mv
from ophyd.signal import Signal


# This is a signal used to indicate if the mirror should be reset at each step
reset_mcm = Signal(name='reset_mcm', value = True)


def reset_mcm_wrapper(plan):
    '''Inserts a 'move mcm to current setpoint' messages before trigger msgs

    This wrapper inserts a msg to move the mcm to it's current setpoints prior
    to each trigger group.
    '''
    trigger_groups = []
    for msg in plan:
        if (msg.command is 'trigger' and
                msg.kwargs['group'] not in trigger_groups):
            trigger_groups.append(msg.kwargs['group'])
            mv_list=[]
            for axis in [mcm.x, mcm.y, mcm.z, mcm.theta, mcm.phi, mcm.chi]:
                mv_list.extend([axis, axis.setpoint.position])
            yield from mv(*mv_list)
        yield msg

reset_mcm_decorator=make_decorator(reset_mcm_wrapper)
