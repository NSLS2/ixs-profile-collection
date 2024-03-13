from ophyd import (Device, Component as Cpt, EpicsMotor, PVPositioner,
                   FormattedComponent as FCpt)

class TempSensor(PVPositioner):
    setpoint = Cpt(EpicsSignal, 'Tset')


class DxtalTemp(Device):
    "Temperatures of the D crystal inside the Analyzer"
    Dtemp = Cpt(TempSensor, 'PID')


D1temp = DxtalTemp('XF:10ID-CT{FbPid:01}', name = 'DxtalTemp')
D2temp = DxtalTemp('XF:10ID-CT{FbPid:02}', name = 'DxtalTemp')
D3temp = DxtalTemp('XF:10ID-CT{FbPid:03}', name = 'DxtalTemp')
D4temp = DxtalTemp('XF:10ID-CT{FbPid:04}', name = 'DxtalTemp')
D5temp = DxtalTemp('XF:10ID-CT{FbPid:05}', name = 'DxtalTemp')
D6temp = DxtalTemp('XF:10ID-CT{FbPid:06}', name = 'DxtalTemp')


# Dtemp1 = EpicsSignal("XF:10ID-CT{FbPid:01}PID.VAL", name="Dtemp1")
# Dtemp2 = EpicsSignal("XF:10ID-CT{FbPid:02}PID.VAL", name="Dtemp2")
# Dtemp3 = EpicsSignal("XF:10ID-CT{FbPid:03}PID.VAL", name="Dtemp3")
# Dtemp4 = EpicsSignal("XF:10ID-CT{FbPid:04}PID.VAL", name="Dtemp4")
# Dtemp5 = EpicsSignal("XF:10ID-CT{FbPid:05}PID.VAL", name="Dtemp5")
# Dtemp6 = EpicsSignal("XF:10ID-CT{FbPid:06}PID.VAL", name="Dtemp6")
