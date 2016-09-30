from ophyd import (EpicsSignalRO, EpicsScaler, DetectorBase, SingleTrigger,
                   Component as Cpt, Device, EpicsSignal, DeviceStatus,
                   ADComponent as ADCpt, Staged)
from ophyd.areadetector.cam import AreaDetectorCam
from ophyd.areadetector import (EpicsSignalWithRBV, ImagePlugin, StatsPlugin,
                                SingleTrigger)


class QuadEM(SingleTrigger, DetectorBase):
    model = Cpt(EpicsSignalRO, 'Model')
    firmware = Cpt(EpicsSignalRO, 'Firmware')

    acquire_mode = Cpt(EpicsSignalWithRBV, 'AcquireMode')
    acquire = Cpt(EpicsSignal, 'Acquire')

    read_format = Cpt(EpicsSignalWithRBV, 'ReadFormat')
    em_range = Cpt(EpicsSignalWithRBV, 'Range')
    ping_pong = Cpt(EpicsSignalWithRBV, 'PingPong')

    integration_time = Cpt(EpicsSignalWithRBV, 'IntegrationTime')
    num_channels = Cpt(EpicsSignalWithRBV, 'NumChannels')
    geometry = Cpt(EpicsSignalWithRBV, 'Geometry')
    resolution = Cpt(EpicsSignalWithRBV, 'Resolution')

    bias_state = Cpt(EpicsSignalWithRBV, 'BiasState')
    bias_interlock = Cpt(EpicsSignalWithRBV, 'BiasInterlock')
    bias_voltage = Cpt(EpicsSignalWithRBV, 'BiasVoltage')
    hvs_readback = Cpt(EpicsSignalRO, 'HVSReadback')
    hvv_readback = Cpt(EpicsSignalRO, 'HVVReadback')
    hvi_readback = Cpt(EpicsSignalRO, 'HVIReadback')

    values_per_read = Cpt(EpicsSignalWithRBV, 'ValuesPerRead')
    sample_time = Cpt(EpicsSignalRO, 'SampleTime_RBV') # yay for consistency
    averaging_time = Cpt(EpicsSignalWithRBV, 'AveragingTime')
    num_average = Cpt(EpicsSignalRO, 'NumAverage_RBV')
    num_averaged = Cpt(EpicsSignalRO, 'NumAveraged_RBV')
    num_acquire = Cpt(EpicsSignalWithRBV, 'NumAcquire')
    num_acquired = Cpt(EpicsSignalRO, 'NumAcquired')
    read_data = Cpt(EpicsSignalRO, 'ReadData')
    ring_overflows = Cpt(EpicsSignalRO, 'RingOverflows')
    trigger_mode = Cpt(EpicsSignal, 'TriggerMode')
    reset = Cpt(EpicsSignal, 'Reset')

    current_name1 = Cpt(EpicsSignal, 'CurrentName1')
    current_name2 = Cpt(EpicsSignal, 'CurrentName2')
    current_name3 = Cpt(EpicsSignal, 'CurrentName3')
    current_name4 = Cpt(EpicsSignal, 'CurrentName4')

    current_offset1 = Cpt(EpicsSignal, 'CurrentOffset1')
    current_offset2 = Cpt(EpicsSignal, 'CurrentOffset2')
    current_offset3 = Cpt(EpicsSignal, 'CurrentOffset3')
    current_offset4 = Cpt(EpicsSignal, 'CurrentOffset4')

    current_offset_calc1 = Cpt(EpicsSignal, 'ComputeCurrentOffset1')
    current_offset_calc2 = Cpt(EpicsSignal, 'ComputeCurrentOffset2')
    current_offset_calc3 = Cpt(EpicsSignal, 'ComputeCurrentOffset3')
    current_offset_calc4 = Cpt(EpicsSignal, 'ComputeCurrentOffset4')

    current_scale1 = Cpt(EpicsSignal, 'CurrentScale1')
    current_scale2 = Cpt(EpicsSignal, 'CurrentScale2')
    current_scale3 = Cpt(EpicsSignal, 'CurrentScale3')
    current_scale4 = Cpt(EpicsSignal, 'CurrentScale4')

    position_offset_x = Cpt(EpicsSignal, 'PositionOffsetX')
    position_offset_y = Cpt(EpicsSignal, 'PositionOffsetY')

    position_offset_calc_x = Cpt(EpicsSignal, 'ComputePosOffsetX')
    position_offset_calc_y = Cpt(EpicsSignal, 'ComputePosOffsetY')

    position_scale_x = Cpt(EpicsSignal, 'PositionScaleX')
    position_scale_Y = Cpt(EpicsSignal, 'PositionScaleY')


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.stage_sigs.update([(self.acquire, 0), # if acquiring, stop
                                (self.acquire_mode, 2) # single mode
                               ])
        self._acquisition_signal = self.acquire


class AH501(QuadEM):
    image = ADCpt(ImagePlugin, 'image1:')
    current1 = ADCpt(StatsPlugin, 'Current1:')
    current2 = ADCpt(StatsPlugin, 'Current2:')
    current3 = ADCpt(StatsPlugin, 'Current3:')
    current4 = ADCpt(StatsPlugin, 'Current4:')


det1 = AH501('XF10ID-BI:AH171:', name='det1')
det2 = AH501('XF10ID-BI:AH172:', name='det2')
det3 = AH501('XF10ID-BI:AH173:', name='det3')
det4 = AH501('XF10ID-BI:AH174:', name='det4')
det5 = AH501('XF10ID-BI:AH175:', name='det5')

for det in [det1, det2, det3, det4, det5]:
    det.configuration_attrs = ['integration_time', 'averaging_time']
    det.read_attrs = ['current1.mean_value','current2.mean_value',
                        'current3.mean_value','current4.mean_value']

sclr = EpicsScaler('XF:10IDD-ES{Sclr:1}', name='sclr')
