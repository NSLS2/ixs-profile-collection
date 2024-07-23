from hkl import E4CV, SimMixin, Lattice, DiffractometerConfiguration
from ophyd import SoftPositioner
from ophyd import Component as Cpt
from hkl.user import *
import pathlib, pyRestTable

class FourCircle(SimMixin, E4CV):
    """
    Our 4-circle.  Eulerian, vertical scattering orientation.
    """
    # the reciprocal axes are defined by SimMixin

#    the = Cpt(SoftPositioner, kind="hinted", init_pos=0)
#    chi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
#    phi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
#    tth = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    th = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:Th}Mtr', labels=('ixs4c',))
    chi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:ChiA}Mtr', labels=('ixs4c',))
    phi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:PhiA}Mtr', labels=('ixs4c',))
    tth = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:2Th}Mtr', labels=('ixs4c',))

def save_config(file_name:str):
# saves the ixs4c_config to a file
    fname = file_name + '.json'
    config_file = config_path / fname
    ixs4c_config.export(config_file)
    print(f"Configuration saved to {config_file}")

def load_config(file_name:str):
# loads the ixs4c_config from a file
    fname = file_name + '.json'
    config_file = config_path / fname
    ixs4c_config.restore(config_file)
    print(f"Configuration loaded from {config_file}")

def preview_config(file_name:str):
# previews the ixs4c_config from a file
    fname = file_name + '.json'
    config_file = config_path / fname
    print(ixs4c_config.preview(config_file))

def all_forward_solutions(hkl_position):
# prints out all forward solutions for the HKLPY
    axes = ixs4c.calc.physical_axis_names
    table = pyRestTable.Table()
    table.labels = axes
    for sol in ixs4c.calc.forward(hkl_position):
        table.addRow([round(getattr(sol, k), 2) for k in axes])

    print(f"solutions for forward({hkl_position}):")
    print(table)

ixs4c = FourCircle("", name="ixs4c")
ixs4c.calc.energy=9.1317
ixs4c.calc.physical_axis_names = {'omega': 'th', 'chi': 'chi', 'phi': 'phi', 'tth': 'tth'}

ixs4c.engine.mode = "constant_phi"
ixs4c.energy.put(9.1317)

ixs4c_config = DiffractometerConfiguration(ixs4c)

config_path = pathlib.Path("/IXS2/data/Ixs4c_config")

## add a sample:
#ixs4c.calc.new_sample('MnF2', lattice=Lattice(a=4.873, b=4.873, c=3.13, alpha=90, beta=90, gamma=90))

## add two, known reflections:
# r1 = ixs4c.calc.sample.add_reflection(1,1,0,position=ixs4c.calc.Position(tth=22.7251, the=7.16641,chi=-4.43721,phi=0))
# r2 = ixs4c.calc.sample.add_reflection(2,1,0,position=ixs4c.calc.Position(tth=36.2993, the=31.8763,chi=1.2541,phi=0))

# add reflections from the current position:
# r3 = ixs4c.calc.sample.add_reflection(1,1,1, ixs4c.real_position}
# 

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

# config_file = config_path / "test_conf.json" -> sets the file name to save the configuration
# ixs4c_config.export(config_file) -> saves the ixs4c configuration
# print(ixs4c_config.preview(config_file)) -> preview the configuration
# ixs4c_config.restore(config_file) -> restores the configuration from the file

# apply constraints:
# e4cv.apply_constraints(
#    {
#       "the": Constraint(10, 40, 0, True),
#       "chi": Constraint(-100, 100, 0, True),
#       "phi": Constraint(-100, 100, 0, True),
#       "tth": Constraint(10, 92.4, 0, True),
#    }
# )

# available modes:
# - bissector
# - constant_chi
# - constant_omega
# - constant_phi
# - double_diffraction
# - psi_constant