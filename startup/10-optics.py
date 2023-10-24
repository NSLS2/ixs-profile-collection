from ophyd import (Device, Component as Cpt, EpicsMotor, PVPositioner,
                   FormattedComponent as FCpt)


# List of available EpicsMotor labels in this script
# [dcm. hrm2, vfm, hfm, xymotor, ssa, table, pinhole, bpm1, bpm1_diag, bpm2, bpm2_diag]


# defined here, also used in 10-optics.py
# why do none of the 'slits' constructs feature gap, center motions?
class Blades(Device):
    top = Cpt(EpicsMotor, '-Ax:T}Mtr')
    bottom = Cpt(EpicsMotor, '-Ax:B}Mtr')
    outboard = Cpt(EpicsMotor, '-Ax:O}Mtr')
    inboard = Cpt(EpicsMotor, '-Ax:I}Mtr')


class DCM(Device):
    y =  Cpt(EpicsMotor, '-Ax:Y}Mtr', labels=('dcm',))
    th = Cpt(EpicsMotor, '-Ax:P}Mtr', labels=('dcm',))
    z2 = Cpt(EpicsMotor, '-Ax:Z2}Mtr', labels=('dcm',))
    p1 = Cpt(EpicsMotor, '-Ax:P1}Mtr', labels=('dcm',))
    r2 = Cpt(EpicsMotor, '-Ax:R2}Mtr', labels=('dcm',))
    pf = Cpt(EpicsMotor, '-Ax:PF}Mtr', labels=('dcm',))


class HRM2(Device):
    ux =  Cpt(EpicsMotor, '-Ax:UX}Mtr', labels=('hrm2',))
    uy =  Cpt(EpicsMotor, '-Ax:UY}Mtr', labels=('hrm2',))
    uth = Cpt(EpicsMotor, '-Ax:UTc}Mtr', labels=('hrm2',))
    uch = Cpt(EpicsMotor, '-Ax:UC}Mtr', labels=('hrm2',))
    uif = Cpt(EpicsMotor, '-Ax:UTI}Mtr', labels=('hrm2',))
    uof = Cpt(EpicsMotor, '-Ax:UTO}Mtr', labels=('hrm2',))
    dx =  Cpt(EpicsMotor, '-Ax:DX}Mtr', labels=('hrm2',))
    dy =  Cpt(EpicsMotor, '-Ax:DY}Mtr', labels=('hrm2',))
    dth = Cpt(EpicsMotor, '-Ax:DTc}Mtr', labels=('hrm2',))
    dch = Cpt(EpicsMotor, '-Ax:DC}Mtr', labels=('hrm2',))
    dif = Cpt(EpicsMotor, '-Ax:DTI}Mtr', labels=('hrm2',))
    dof = Cpt(EpicsMotor, '-Ax:DTO}Mtr', labels=('hrm2',))
    d1 =  Cpt(EpicsMotor, '-Pico:m1}Mtr', labels=('hrm2',))
    d2 =  Cpt(EpicsMotor, '-Pico:m2}Mtr', labels=('hrm2',))
    d3 =  Cpt(EpicsMotor, '-Pico:m3}Mtr', labels=('hrm2',))
    d4 =  Cpt(EpicsMotor, '-Pico:m4}Mtr', labels=('hrm2',))
    d5 =  Cpt(EpicsMotor, '-Pico:m5}Mtr', labels=('hrm2',))
    bs =  Cpt(EpicsMotor, '-Pico:m6}Mtr', labels=('hrm2',))



class VFM(Device):
    ux = Cpt(EpicsMotor, '-Ax:UX}Mtr', labels=('vfm',))
    uy = Cpt(EpicsMotor, '-Ax:UY}Mtr', labels=('vfm',))
    dx = Cpt(EpicsMotor, '-Ax:DX}Mtr', labels=('vfm',))
    dy = Cpt(EpicsMotor, '-Ax:DY}Mtr', labels=('vfm',))
    ub = Cpt(EpicsMotor, '-Ax:UB}Mtr', labels=('vfm',))
    db = Cpt(EpicsMotor, '-Ax:DB}Mtr', labels=('vfm',))


class HFM(Device):
    ux = Cpt(EpicsMotor, '-Ax:UX}Mtr', labels=('hfm',))
    uy = Cpt(EpicsMotor, '-Ax:UY}Mtr', labels=('hfm',))
    dx = Cpt(EpicsMotor, '-Ax:DX}Mtr', labels=('hfm',))
    dy = Cpt(EpicsMotor, '-Ax:DY}Mtr', labels=('hfm',))
    ub = Cpt(EpicsMotor, '-Ax:UB}Mtr', labels=('hfm',))
    db = Cpt(EpicsMotor, '-Ax:DB}Mtr', labels=('hfm',))


class XYMotor(Device):
    x = Cpt(EpicsMotor, '-Ax:X}Mtr', labels=('xymotor',))
    y = Cpt(EpicsMotor, '-Ax:Y}Mtr', labels=('xymotor',))


class SSA(Device):
    top = Cpt(EpicsMotor, '-Ax:T}Mtr', labels=('ssa',))
    bottom = Cpt(EpicsMotor, '-Ax:B}Mtr', labels=('ssa',))


class Table(Device):
    x = Cpt(EpicsMotor, '-Ax:X4}Mtr', labels=('table',))
    y = Cpt(EpicsMotor, '-Ax:Y4}Mtr', labels=('table',))
    th =Cpt(EpicsMotor, '-Ax:X4a1}Mtr', labels=('table',))


class Pinhole(Device):
    ux = Cpt(EpicsMotor, '-Ax:UX}Mtr', labels=('pinhole',))
    uy = Cpt(EpicsMotor, '-Ax:UY}Mtr', labels=('pinhole',))
    dx = Cpt(EpicsMotor, '-Ax:DX}Mtr', labels=('pinhole',))
    dy = Cpt(EpicsMotor, '-Ax:DY}Mtr', labels=('pinhole',))


class MCMBase(PVPositioner):
    setpoint = FCpt(EpicsSignal, '{self.prefix}{self._ch_name}')
    readback = FCpt(EpicsSignal, '{self.prefix}{self._ch_name}')
    actuate = Cpt(EpicsSignal, '}Mov')
    actuate_value = 1
    stop_signal = Cpt(EpicsSignal, '}Kill')
    stop_value = 1
    # all six axes are coupled, so 'InPos' is an array of six values
    done = Cpt(EpicsSignal, '}InPos')
    done_value = True

    def __init__(self, prefix, ch_name=None, **kwargs):
        self._ch_name = ch_name
        super().__init__(prefix, **kwargs)

    @property
    def moving(self):
        return (self.done_value != all(self.done.get(use_monitor=False)))


class MCM(Device):
    x = Cpt(MCMBase, '', ch_name='-Ax:X}Mtr')
    y = Cpt(MCMBase, '', ch_name='-Ax:Y}Mtr')
    z = Cpt(MCMBase, '', ch_name='-Ax:Z}Mtr')
    theta = Cpt(MCMBase, '', ch_name='-Ax:Rx}Mtr')
    phi = Cpt(MCMBase, '', ch_name='-Ax:Ry}Mtr')
    chi = Cpt(MCMBase, '', ch_name='-Ax:Rz}Mtr')


dcm = DCM('XF:10IDA-OP{Mono:DCM', name='dcm')
vfm = VFM('XF:10IDD-OP{VFM:1', name='vfm')
hfm = HFM('XF:10IDD-OP{HFM:1', name='hfm')
hrm2 = HRM2('XF:10IDB-OP{Mono:HRM2', name='hrm2')

s1 = Blades('XF:10IDA-OP{Slt:1', name='s1')
s2 = Blades('XF:10IDC-OP{Slt:4', name='s2')
s3 = Blades('XF:10IDD-OP{Slt:5', name='s3')

bpm1 = XYMotor('XF:10IDA-OP{BPM:1', name='bpm1', labels=('bpm1',))
bpm1_diag = EpicsMotor('XF:10IDA-BI{BPM:1-Ax:YFoil}Mtr', name='bpm1_diag', labels=('bpm1_diag',))
bpm2 = XYMotor('XF:10IDC-OP{BPM:2', name='bpm2', labels=('bpm2',))
bpm2_diag = EpicsMotor('XF:10IDC-BI{BPM:2-Ax:Y}Mtr', name='bpm2_diag', labels=('bpm2_diag',))

ssa = SSA('XF:10IDB-OP{SSA:1', name='ssa')
k3 = Table('XF:10IDC-OP{Tbl:1', name='k3')
ph = Pinhole('XF:10IDD-OP{Pinh:1', name='ph')

mcm = MCM('XF:10IDD-OP{MCM:1', name='mcm')
