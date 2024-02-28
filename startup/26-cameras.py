from ophyd import ProsilicaDetector, SingleTrigger

from ophyd import Component as Cpt

from nslsii.ad33 import StatsPluginV33 as StatsPlugin


class StandardProsilica(SingleTrigger, ProsilicaDetector):
    stats1 = Cpt(StatsPlugin, 'Stats1:')
    stats2 = Cpt(StatsPlugin, 'Stats2:')
    stats3 = Cpt(StatsPlugin, 'Stats3:')
    stats4 = Cpt(StatsPlugin, 'Stats4:')
    stats5 = Cpt(StatsPlugin, 'Stats5:')
    

cam1 = StandardProsilica('XF:10IDA-BI{BPM1-Cam:1}', name='cam1')
cam2 = StandardProsilica('XF:10IDC-BI{BPM2-Cam:1}', name='cam2')