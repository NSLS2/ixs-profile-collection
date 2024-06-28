from hkl import E4CV, SimMixin, Lattice, DiffractometerConfiguration
from ophyd import SoftPositioner
from ophyd import Component as Cpt
from hkl.user import *

class FourCircle(SimMixin, E4CV):
    """
    Our 4-circle.  Eulerian, vertical scattering orientation.
    """
    # the reciprocal axes are defined by SimMixin

    the = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    chi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    phi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    tth = Cpt(SoftPositioner, kind="hinted", init_pos=0)

ixs4c = FourCircle("", name="ixs4c")
ixs4c.calc.physical_axis_names = {'omega': 'the', 'chi': 'chi', 'phi': 'phi', 'tth': 'tth'}
ixs4c.engine.mode = "constant_phi"
ixs4c.energy.put(9.1317)

## aa a sample
#ixs4c.calc.new_sample('MnF2', lattice=Lattice(a=4.873, b=4.873, c=3.13, alpha=90, beta=90, gamma=90))

## add two, known reflections
# r1 = ixs4c.calc.sample.add_reflection(1,1,0,position=ixs4c.calc.Position(tth=22.7251, the=7.16641,chi=-4.43721,phi=0))
# r2=ixs4c.calc.sample.add_reflection(2,1,0,position=ixs4c.calc.Position(tth=36.2993, the=31.8763,chi=1.2541,phi=0))

## compute orientation matrix UB
# ixs4c.calc.sample.compute_UB(r1,r2)

## display the setup
# ixs4c.pa()

## forward calculation of motor positions from HKL
# ixs4c.forward((1, 0.5, 0))

## inverse calculation of HKL from motor positions
# ixs4c.inverse((the_pos, chi_pos, phi_pos, tth_pos))

## Useful basic commands
# ixs4c.position -> HKL related to the current physical position
# ixs4c.real_position -> current physical position
# ixs4c.calc.forward((1,0,1)) -> calculates physical motor positions
# given an HKL
# ixs4c.move(0,0,0) -> physical motion corresponding to a given HKL