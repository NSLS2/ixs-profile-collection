from ophyd import Component as Cpt, Device, EpicsSignal


class URATemperature(Device):
    d1temp = Cpt(EpicsSignal, '01}PID.VAL', labels=("DxtalTemp", ))
    d2temp = Cpt(EpicsSignal, '02}PID.VAL', labels=("DxtalTemp", ))
    d3temp = Cpt(EpicsSignal, '03}PID.VAL', labels=("DxtalTemp", ))
    d4temp = Cpt(EpicsSignal, '04}PID.VAL', labels=("DxtalTemp", ))
    d5temp = Cpt(EpicsSignal, '05}PID.VAL', labels=("DxtalTemp", ))
    d6temp = Cpt(EpicsSignal, '06}PID.VAL', labels=("DxtalTemp", ))

ura_temp = URATemperature("XF:10IDD-CT{EPid:", name="uratemperature")