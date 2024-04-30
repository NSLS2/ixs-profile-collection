from ophyd import (Component as Cpt, Device, EpicsMotor, EpicsSignal)

# List of available EpicsMotor labels in this script
# [analyzer, spectrometer, analyzerdxtals, analyzerslits, mcmslits, samplestage, whl, anapd, anpd]

class Analyzer(Device):
    uy = Cpt(EpicsMotor, '-OP{Analy:1-Ax:UY}Mtr', labels=('analyzer',))
    dy = Cpt(EpicsMotor, '-OP{Analy:1-Ax:DY}Mtr', labels=('analyzer',))
    ay = Cpt(EpicsMotor, '-OP{Analy:1-Ax:A}Mtr', labels=('analyzer',))
    by = Cpt(EpicsMotor, '-OP{Analy:1-Ax:B}Mtr', labels=('analyzer',))

    cfth = Cpt(EpicsMotor, '-ES{Ez4:1-Ax:1}Mtr', labels=('analyzer',))
    cchi = Cpt(EpicsMotor, '-ES{Ez4:1-Ax:2}Mtr', labels=('analyzer',))
    wfth = Cpt(EpicsMotor, '-ES{Ez4:1-Ax:3}Mtr', labels=('analyzer',))
    wchi = Cpt(EpicsMotor, '-ES{Ez4:1-Ax:4}Mtr', labels=('analyzer',))


class Spectrometer(Device):
    tth = Cpt(EpicsMotor, '-OP{Spec:1-Ax:2Th}Mtr', labels=('spectrometer',))
    th = Cpt(EpicsMotor, '-OP{Spec:1-Ax:Th}Mtr', labels=('spectrometer',))
    chi = Cpt(EpicsMotor, '-OP{Spec:1-Ax:ChiA}Mtr', labels=('spectrometer',))
    phi = Cpt(EpicsMotor, '-OP{Spec:1-Ax:PhiA}Mtr', labels=('spectrometer',))


class AnalyzerDXtals(Device):
    d1the = Cpt(EpicsMotor, '2-Ax:1}Mtr', labels=('analyzerdxtals',))
    d1phi = Cpt(EpicsMotor, '2-Ax:2}Mtr', labels=('analyzerdxtals',))
    d2the = Cpt(EpicsMotor, '2-Ax:3}Mtr', labels=('analyzerdxtals',))
    d2phi = Cpt(EpicsMotor, '2-Ax:4}Mtr', labels=('analyzerdxtals',))
    d3the = Cpt(EpicsMotor, '3-Ax:1}Mtr', labels=('analyzerdxtals',))
    d3phi = Cpt(EpicsMotor, '3-Ax:2}Mtr', labels=('analyzerdxtals',))
    d4the = Cpt(EpicsMotor, '3-Ax:3}Mtr', labels=('analyzerdxtals',))
    d4phi = Cpt(EpicsMotor, '3-Ax:4}Mtr', labels=('analyzerdxtals',))
    d5the = Cpt(EpicsMotor, '4-Ax:1}Mtr', labels=('analyzerdxtals',))
    d5phi = Cpt(EpicsMotor, '4-Ax:2}Mtr', labels=('analyzerdxtals',))
    d6the = Cpt(EpicsMotor, '4-Ax:3}Mtr', labels=('analyzerdxtals',))
    d6phi = Cpt(EpicsMotor, '4-Ax:4}Mtr', labels=('analyzerdxtals',))
    anpd = Cpt(EpicsMotor,  '5-Ax:1}Mtr', labels=('analyzerdxtals',))


class AnalyzerSlits(Device):
    top = Cpt(EpicsMotor,  '5-Ax:2}Mtr', labels=('analyzerslits',))
    bottom = Cpt(EpicsMotor,  '5-Ax:3}Mtr', labels=('analyzerslits',))
    outboard = Cpt(EpicsMotor,  '7-Ax:3}Mtr', labels=('analyzerslits',))
    inboard = Cpt(EpicsMotor,  '7-Ax:4}Mtr', labels=('analyzerslits',))


class MCMSlits(Device):
 #   top = Cpt(EpicsMotor, '6-Ax:3}Mtr', labels=('mcmslits',))
 #   bottom = Cpt(EpicsMotor, '6-Ax:4}Mtr', labels=('mcmslits',))
    inboard = Cpt(EpicsMotor,  '-Ax:Xi}Mtr', labels=('mcmslits',))
    outboard = Cpt(EpicsMotor,  '-Ax:Xo}Mtr', labels=('mcmslits',))


class SampleStage(Device):
    ty = Cpt(EpicsMotor, '{Spec:1-Ax:Y}Mtr', labels=('samplestage',))
    tx = Cpt(EpicsMotor, '{Spec:1-Ax:X}Mtr', labels=('samplestage',))
    tz = Cpt(EpicsMotor, '{Spec:1-Ax:Z}Mtr', labels=('samplestage',))
    sy = Cpt(EpicsMotor, '{Env:1-Ax:Y}Mtr', labels=('samplestage',))
    sx = Cpt(EpicsMotor, '{Env:1-Ax:X}Mtr', labels=('samplestage',))
    sz = Cpt(EpicsMotor, '{Env:1-Ax:Z}Mtr', labels=('samplestage',))


analyzer = Analyzer('XF:10IDD', name='analyzer')
spec = Spectrometer('XF:10IDD', name='spec')
analyzer_xtals = AnalyzerDXtals('XF:10IDD-ES{Ez4:', name='analyzer_xtals')
analyzer_slits = AnalyzerSlits('XF:10IDD-ES{Ez4:', name='analyzer_slits')
mcm_slits = MCMSlits('XF:10IDD-OP{MCMSlt1', name='mcmslits')
sample_stage = SampleStage('XF:10IDD-OP', name='s')

whl = EpicsMotor('XF:10IDD-OP{Abs:1-Ax:Wheel}Mtr', name='whl', labels=('whl',))

anapd = EpicsMotor('XF:10IDD-ES{Ez4:8-Ax:3}Mtr', name='anapd', labels=('anapd',))
anpd = EpicsMotor('XF:10IDD-ES{Ez4:5-Ax:1}Mtr', name='anpd', labels=('anpd',))
airpad = EpicsSignal("XF:10IDD-CT{IOC-MC:12}AirOn-cmd", name="airpad")

# name_you_want = EpicsMotor(PVNAME_BASE, name='name_you_want')