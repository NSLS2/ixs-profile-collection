from ophyd.areadetector.cam import Lambda750kCam
from ophyd import Component as Cpt
from ophyd import ROIPlugin, TransformPlugin, EpicsSignal, EpicsSignalRO
from ophyd.areadetector.base import EpicsSignalWithRBV as SignalWithRBV
from ophyd.areadetector.cam import CamBase
from ophyd.areadetector import ADComponent as ADCpt, DetectorBase
from nslsii.ad33 import StatsPluginV33, SingleTriggerV33
from ophyd.areadetector.plugins import PluginBase


class PluginCV(PluginBase):
    ...

class LambdaDetector(DetectorBase):
    _html_docs = ['lambda.html']
    cam = Cpt(Lambda750kCam, 'cam1:')

class Lambda(SingleTriggerV33, LambdaDetector):
    # MR20200122: created all dirs recursively in /nsls2/jpls/data/lambda/
    # from 2020 to 2030 with 777 permissions, owned by xf12id1 user.
    # tiff = Cpt(TIFFPluginWithFileStore,
    #            suffix="TIFF1:",
    #            write_path_template="/nsls2/xf12id1g/data/lambda/%Y/%m/%d/",
    #            read_path_template="/nsls2/xf12id1g/data/lambda/%Y/%m/%d/",
    #            root='/nsls2/xf12id1g/data') 
    cv1 = Cpt(PluginCV, 'CV1:')
    
    roi1 = Cpt(ROIPlugin, 'ROI1:')
    roi2 = Cpt(ROIPlugin, 'ROI2:')
    roi3 = Cpt(ROIPlugin, 'ROI3:')
    roi4 = Cpt(ROIPlugin, 'ROI4:')
    roi5 = Cpt(ROIPlugin, 'ROI5:')
    roi6 = Cpt(ROIPlugin, 'ROI6:')
    roi7 = Cpt(ROIPlugin, 'ROI7:')


    stats1 = Cpt(StatsPluginV33, 'Stats1:', read_attrs=['total'])
    stats2 = Cpt(StatsPluginV33, 'Stats2:', read_attrs=['total'])
    stats3 = Cpt(StatsPluginV33, 'Stats3:', read_attrs=['total'])
    stats4 = Cpt(StatsPluginV33, 'Stats4:', read_attrs=['total'])
    stats5 = Cpt(StatsPluginV33, 'Stats5:', read_attrs=['total'])
    stats6 = Cpt(StatsPluginV33, 'Stats6:', read_attrs=['total'])
    stats7 = Cpt(StatsPluginV33, 'Stats7:', read_attrs=['total'])


    trans1 = Cpt(TransformPlugin, 'Trans1:')

    low_thr = Cpt(EpicsSignal, 'cam1:LowEnergyThreshold')
    hig_thr = Cpt(EpicsSignal, 'cam1:HighEnergyThreshold')
    oper_mode = Cpt(EpicsSignal, 'cam1:OperatingMode')

lambda_det = Lambda('XF:10IDC-BI{Lambda-Cam:1}', name='lambda_det')
for j in range(1, 8):
    getattr(lambda_det, f'stats{j}').kind = 'normal'
lambda_det.stats7.total.kind = 'hinted'


# Impose Stats4 to be ROI4 if in the future we need to exclude bad pixels
def set_defaut_stat_roi():
    yield from bps.mv(lambda_det.stats1.nd_array_port, 'ROI1')
    yield from bps.mv(lambda_det.stats2.nd_array_port, 'ROI2')
    yield from bps.mv(lambda_det.stats3.nd_array_port, 'ROI3')
    yield from bps.mv(lambda_det.stats4.nd_array_port, 'ROI4')


def set_lambda_exposure(exposure):
    # Sets the Lambda detector exposure time (exposure)
    det = lambda_det
    yield from bps.mv(det.cam.acquire_time, exposure, det.cam.acquire_period, exposure)

