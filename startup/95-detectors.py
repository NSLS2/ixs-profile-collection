from ophyd import (EpicsSignalRO, EpicsScaler, DetectorBase, SingleTrigger,
                   Component as Cpt, Device, EpicsSignal, DeviceStatus,
                   ADComponent as ADCpt, Staged)
from ophyd.areadetector.cam import AreaDetectorCam
from ophyd.areadetector import (EpicsSignalWithRBV, ImagePlugin, StatsPlugin,
                                SingleTrigger)


from ophyd.quadem import NSLS_EM, TetrAMM, QuadEMPort

det1 = NSLS_EM('XF10ID-BI:AH171:', name='det1')
det2 = NSLS_EM('XF10ID-BI:AH172:', name='det2')
det3 = NSLS_EM('XF10ID-BI:AH173:', name='det3')
det4 = NSLS_EM('XF10ID-BI:AH174:', name='det4')
det5 = NSLS_EM('XF10ID-BI:AH175:', name='det5')

for j, det in enumerate([det1, det2, det3, det4, det5]):
    det.configuration_attrs = ['integration_time', 'averaging_time']
    det.read_attrs = ['current1.mean_value','current2.mean_value',
                        'current3.mean_value','current4.mean_value']
    det.conf.port_name.put(f'AH17{j+1}')

tm1 = NSLS_EM('XF:10ID-BI:TM176:', name='tm1')
tm2 = NSLS_EM('XF:10ID-BI:TM178:', name='tm2')

for j, det in enumerate([tm1, tm2]):
    det.configuration_attrs = ['integration_time', 'averaging_time']
    det.read_attrs = ['current1.mean_value','current2.mean_value',
                        'current3.mean_value','current4.mean_value']
    det.conf.port_name.put(f'TM17{2*j+6}')

sclr = EpicsScaler('XF:10IDD-ES{Sclr:1}', name='sclr')
for j in range(1, 33):
    getattr(sclr.channels, f'chan{j}').kind= 'omitted'

sclr.channels.chan21.name = 'I0'
sclr.channels.chan21.kind = 'hinted'


    