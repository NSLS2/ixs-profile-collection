'''
from bluesky.callbacks.olog import OlogCallback
from bluesky.callbacks import LivePlot


from suitcase.spec import DocumentToSpec

spec_cb = DocumentToSpec('/tmp/spec1.spec')
RE.subscribe(spec_cb)
'''

from utils.sixcircle import SixCircle

sc = SixCircle()

def whsc():
# Defines the current position of the sample
    Th = spec.th.position
    Tth = spec.tth.position
    Chi = spec.chi.position
    Phi = spec.phi.position
    sc.mv(tth=Tth,th=Th,chi=Chi,phi=Phi)


def wh():
    whsc()
    # update_spec_signals()


def wh_refresh():
    whsc()
    return sc.H, sc.K, sc.L


def freeze(*arg):
    sc.freeze(*arg)


def setfrozen(*arg):
    sc.setfrozen(*arg)


def pa():
    sc.pa()


def or0(*args):
# defines the OR0 vector from the current angle positions
    whsc()
    sc.or0(*args)


def or1(*args):
# defines the OR0 vector from the current angle positions
    whsc()
    sc.or1(*args)

# def update_spec_signals():
#    spec.H.put(sc.H)
#    spec.K.put(sc.K)
#    spec.L.put(sc.L)
#    spec.Q.put(sc.ABSQ*10)
#    spec.LAMBDA.put(sc.LAMBDA)
#    spec.HAZ.put(sc.g_haz)
#    spec.KAZ.put(sc.g_kaz)
#    spec.LAZ.put(sc.g_laz)
#    spec.AZIMUTH.put(sc.AZIMUTH)
#    spec.ALPHA.put(sc.ALPHA)
#    spec.BETA.put(sc.BETA)
#    spec.OMEGA.put(sc.OMEGA)


def br(*args):
    sc.br(*args)
    # update_spec_signals()
    yield from bps.mv(spec.th, sc.TH, spec.tth, sc.TTH, spec.chi, sc.CHI, spec.phi, sc.PHI)

# def br(h, k, l):
#     sc.br(h, k, l)
#     RE.md.update({'H', h, 'K', k, 'L', l})
#     yield from bps.mv(spec.th, sc.TH, spec.tth, sc.TTH, spec.chi, sc.CHI, spec.phi, sc.PHI)

def sc_init(conf_file):
    sc.load(conf_file)
    setfrozen(345)
    freeze(0,-0.401,0)
    sc.mv(mu=-0.401)
    wh()
#    update_spec_signals()


def ca(h, k, l):
    res = sc.ca(h,k,l)


#def mv()
#    if mov == True and result is not None:
#        sc.mv(tth = result[0], th = result[1], chi = result[2], phi = result[3])


sc_init('/IXS2/data/run2025/dia_test.conf')
