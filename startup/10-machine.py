from ophyd import (Device, EpicsSignal, EpicsSignalRO,
                   EpicsMotor, PVPositioner, PVPositionerPC,
                   Component as Cpt)


# List of available EpicsMotor labels in this script
# [crl. fes]

# SR current
# CNT012 src  SRcur  SR:C03-BI{DCCT:1}I:Real-I
sr_curr = EpicsSignalRO('SR:OPS-BI{DCCT:1}I:Real-I', name='sr_curr')

class CRL(Device):
    x = Cpt(EpicsMotor, '-Ax:X}Mtr', labels=('crl',))
    y = Cpt(EpicsMotor, '-Ax:Y}Mtr', labels=('crl',))
    th =Cpt(EpicsMotor, '-Ax:P}Mtr', labels=('crl',))


class Slit2DGap(PVPositionerPC):
    readback = Cpt(EpicsSignalRO, 't2.C')
    setpoint = Cpt(EpicsSignal, 'size')
    done = Cpt(EpicsSignalRO, 'DMOV')
    done_value = 1


class Slit2DCenter(PVPositionerPC):
    readback = Cpt(EpicsSignalRO, 't2.D')
    setpoint = Cpt(EpicsSignal, 'center')
    done = Cpt(EpicsSignalRO, 'DMOV')
    done_value = 1


class Slit2D(Device):
    "Center and gap with virtual motors"
    xc = Cpt(Slit2DCenter, '12-Ax:X}')
    yc = Cpt(Slit2DCenter, '12-Ax:Y}')
    xg = Cpt(Slit2DGap, '12-Ax:X}')
    yg = Cpt(Slit2DGap, '12-Ax:Y}')


class Slit2DBlades(Device):
    top = Cpt(EpicsMotor, '1-Ax:T}Mtr', labels=('fes',))
    bottom = Cpt(EpicsMotor, '2-Ax:B}Mtr', labels=('fes',))
    outboard = Cpt(EpicsMotor, '1-Ax:O}Mtr', labels=('fes',))
    inboard = Cpt(EpicsMotor, '2-Ax:I}Mtr', labels=('fes',))


class FESlits(Slit2DBlades, Slit2D):
    "combine t b i o and xc yc xg yg"
    pass


# TODO: revisit the IVU as a PseudoPositioner
class Undulator(PVPositioner):
    # setpoint = Cpt(EpicsSignal, '}Man:SP:Gap')
    setpoint = Cpt(EpicsSignal, '-Ax:Gap}-Mtr-SP')
    # Note: there are actually 2 readbacks for gap position, Y1 & Y2.
    # This should be fixed at the EPICS level to provide an avg gap
    # readback = Cpt(EpicsSignalRO, '}Y1:Rbv')
    readback = Cpt(EpicsSignalRO, '-Ax:Gap}-Mtr.RBV')
    # actuate = Cpt(EpicsSignal, '}ManG:Go_.PROC')
    actuate = Cpt(EpicsSignal, '-Ax:Gap}-Mtr-Go')
    # done = Cpt(EpicsSignalRO, '-Mtr:Gap}.DMOV')
    done = Cpt(EpicsSignalRO, '-Ax:Gap}-Mtr.DMOV')
    # stop_signal = Cpt(EpicsSignal, '}Man:Stop_.PROC')
    stop_signal = Cpt(EpicsSignal, '-Ax:Gap}-Mtr.STOP')

    actuate_value = 1
    done_value = 1
    stop_value = 1


class SROFB(Device):
    uofb_pv = Cpt(EpicsSignal, '}ConfigMode-I')
    id_bump_pv = Cpt(EpicsSignal, 'C10-ID}Enabled-I')
    nudge_pv = Cpt(EpicsSignal, 'C10-ID}Nudge-Enabled')
    nudge_increment = Cpt(EpicsSignal, 'C10-ID}angle-increment-SP')
    horz_plane_nudge = Cpt(EpicsSignal, 'C10-ID}Nudge:X')
    vert_plane_nudge = Cpt(EpicsSignal, 'C10-ID}Nudge:Y')
    nudge_status = Cpt(EpicsSignal, 'C10-ID}Nudge-StatusMsg')
    xa_rbv = Cpt(EpicsSignal, 'BUMP:C10-IXS}angle:X-SP')
    xy_rbv = Cpt(EpicsSignal, 'BUMP:C10-IXS}angle:Y-SP')
    

ivu22 = Undulator('SR:C10-ID:G1{IVU22:1', name='ivu22')
ivu22.readback.name = 'ivu22'
fes = FESlits('FE:C10A-OP{Slt:', name='fes')
crl = CRL('FE:C10A-OP{CRL:1', name='crl')

strg_ring_orb_fb = SROFB('SR:UOFB{', name='srofb')
