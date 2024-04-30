from bluesky.suspenders import SuspendFloor
from ophyd import EpicsSignal

tm1sum = EpicsSignal('XF:10ID-BI:TM176:SumAll:MeanValue_RBV')
susp = SuspendFloor(tm1sum, 1.e-5, resume_thresh = 1.e-5, sleep = 1*60)
