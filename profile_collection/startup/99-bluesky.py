from bluesky.callbacks.olog import OlogCallback
from bluesky.callbacks import LivePlot
from bluesky.global_state import gs

from bluesky.callbacks.scientific import plot_peak_stats


# gs.DETS.append(det4)

#from bluesky.scientific_callbacks import plot_peak_stats

from bluesky.plans import *

from bluesky.spec_api import ct, ascan, dscan


from suitcase.spec import DocumentToSpec

spec_cb = DocumentToSpec('/tmp/spec1.spec')

import bluesky.spec_api
from bluesky.plans import planify, subs_context
from functools import wraps

def print_ps_summary(ps=None, plot=False):
    if ps is None:
        ps = gs.PS

    print('cen: {ps.cen}\nmax_pos: {ps.max}\nfwhm: {ps.fwhm}'.format(ps=ps))
    if plot:
        plot_peak_stats(gs.PS)

@wraps(bluesky.spec_api.dscan)
@planify
def dscan(*args, **kwargs):
    plans = []
    with subs_context(plans, {'all': [spec_cb]}):
        plans.append(bluesky.spec_api.dscan(*args, **kwargs))
    return plans


@wraps(bluesky.spec_api.ascan)
@planify
def ascan(*args, **kwargs):
    plans = []
    with subs_context(plans, [spec_cb]):
        plans.append(bluesky.spec_api.ascan(*args, **kwargs))
    return plans


@wraps(bluesky.spec_api.ct)
@planify
def ct(*args, **kwargs):
    plans = []
    with subs_context(plans, [spec_cb]):
        plans.append(bluesky.spec_api.ct(*args, **kwargs))
    return plans


#relabel_motors()
def scan_2_slts_fe(start, stop):
     ret1 = yield from dscan(slt1, start, stop)
     st2, stp2 = some_comp(ret1)
     ret2 = yield from dscan(slt2, st2, stp2)

def count_around(plan):
    yield from ct()
    yield from plan
    yield from ct()
