from hkl import E4CV, SimMixin
from ophyd import SoftPositioner
from ophyd import Component as Cpt

class FourCircle(SimMixin, E4CV):
    """
    Our 4-circle.  Eulerian, vertical scattering orientation.
    """
    # the reciprocal axes are defined by SimMixin

    omega = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    chi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    phi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    tth = Cpt(SoftPositioner, kind="hinted", init_pos=0)

ixs4c = FourCircle("", name="ixs4c")
