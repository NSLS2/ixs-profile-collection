import bluesky.plans as bp
import bluesky.preprocessors as bpp
import bluesky.plan_stubs as bps
from bluesky.callbacks.fitting import PeakStats


def dcm_setup():
    det = tm1
    yname = tm1.sum_all.mean_value.name
    ps = PeakStats(dcm.p1.user_readback.name, yname)
    yield from bpp.subs_wrapper(bp.rel_scan([det], dcm.p1, -80, 80, 40), ps)

    cen = ps.cen
    com = ps.com
    fwhm = ps.fwhm

    print('Peak stats\n', ps)
    if fwhm < 50 and abs(cen - com)/ fwhm < 1 and len(ps.crossings) == 2:
        yield from bps.mv(dcm.p1, cen)
        print("DCM moved to center")
    else:
        print("Do not think we found a peak. Motor not moved!")


def hrm_in():
    yield from bps.mv(hrm2.ux, 0, hrm2.dx, 0, hrm2.bs, 3)


def hrm_out():
    yield from bps.mv(hrm2.ux, -20, hrm2.dx, -20, hrm2.bs, 0)