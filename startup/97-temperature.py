from ophyd import EpicsSignal


d1temp = EpicsSignal("XF:10ID-CT{FbPid:01}PID.VAL", name="d1temp", labels=("DxtalTemp", ))
d2temp = EpicsSignal("XF:10ID-CT{FbPid:02}PID.VAL", name="d2temp", labels=("DxtalTemp", ))
d3temp = EpicsSignal("XF:10ID-CT{FbPid:03}PID.VAL", name="d3temp", labels=("DxtalTemp", ))
d4temp = EpicsSignal("XF:10ID-CT{FbPid:04}PID.VAL", name="d4temp", labels=("DxtalTemp", ))
d5temp = EpicsSignal("XF:10ID-CT{FbPid:05}PID.VAL", name="d5temp", labels=("DxtalTemp", ))
d6temp = EpicsSignal("XF:10ID-CT{FbPid:06}PID.VAL", name="d6temp", labels=("DxtalTemp", ))
