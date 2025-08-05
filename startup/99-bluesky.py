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
    sc.TH = spec.th.position
    sc.TTH = spec.tth.position
    sc.CHI = spec.chi.position
    sc.PHI = spec.phi.position


def wh():
    whsc()
    sc.wh_refresh()
    sc.wh()


def or0(*args):
# defines the OR0 vector from the current angle positions
    whsc()
    sc.or0(*args)


def spec_mv(*args):
    sc.br(*args)
    yield RE(bps.mvr(spec.th, sc.TH, spec.tth, sc.TTH, spec.chi, sc.CHI, spec.phi, sc.PHI))