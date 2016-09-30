from ophyd import (Component as Cpt, Device, EpicsMotor)


class Analyzer(Device):
    uy = Cpt(EpicsMotor, '-OP{Analy:1-Ax:UY}Mtr')
    dy = Cpt(EpicsMotor, '-OP{Analy:1-Ax:DY}Mtr')
    ay = Cpt(EpicsMotor, '-OP{Analy:1-Ax:A}Mtr')
    by = Cpt(EpicsMotor, '-OP{Analy:1-Ax:B}Mtr')
    tth = Cpt(EpicsMotor, '-OP{Spec:1-Ax:2Th}Mtr')
    th =  Cpt(EpicsMotor, '-OP{Spec:1-Ax:Th}Mtr')
    chi = Cpt(EpicsMotor, '-OP{Spec:1-Ax:ChiA}Mtr')
    phi = Cpt(EpicsMotor, '-OP{Spec:1-Ax:PhiA}Mtr')
    cfth = Cpt(EpicsMotor, '-ES{Ez4:1-Ax:1}Mtr')
    cchi = Cpt(EpicsMotor, '-ES{Ez4:1-Ax:2}Mtr')
    wfth = Cpt(EpicsMotor, '-ES{Ez4:1-Ax:3}Mtr')
    wchi = Cpt(EpicsMotor, '-ES{Ez4:1-Ax:4}Mtr')


class AnalyzerDXtals(Device):
    d1the = Cpt(EpicsMotor, '2-Ax:1}Mtr')
    d1phi = Cpt(EpicsMotor, '2-Ax:2}Mtr')
    d2the = Cpt(EpicsMotor, '2-Ax:3}Mtr')
    d2phi = Cpt(EpicsMotor, '2-Ax:4}Mtr')
    d3the = Cpt(EpicsMotor, '3-Ax:1}Mtr')
    d3phi = Cpt(EpicsMotor, '3-Ax:2}Mtr')
    d4the = Cpt(EpicsMotor, '3-Ax:3}Mtr')
    d4phi = Cpt(EpicsMotor, '3-Ax:4}Mtr')
    d5the = Cpt(EpicsMotor, '4-Ax:1}Mtr')
    d5phi = Cpt(EpicsMotor, '4-Ax:2}Mtr')
    d6the = Cpt(EpicsMotor, '4-Ax:3}Mtr')
    d6phi = Cpt(EpicsMotor, '4-Ax:4}Mtr')
    anpd = Cpt(EpicsMotor,  '5-Ax:1}Mtr')


class AnalyzerSlits(Device):
    top = Cpt(EpicsMotor,  '5-Ax:2}Mtr')
    bottom = Cpt(EpicsMotor,  '5-Ax:3}Mtr')

class MCMSlits(Device):
    top = Cpt(EpicsMotor, '6-Ax:3}Mtr')
    bottom = Cpt(EpicsMotor, '6-Ax:4}Mtr')
    inboard = Cpt(EpicsMotor,  '7-Ax:1}Mtr')
    outboard = Cpt(EpicsMotor,  '7-Ax:2}Mtr')


class SampleStage(Device):
    ty = Cpt(EpicsMotor, '{Spec:1-Ax:Y}Mtr')
    tx = Cpt(EpicsMotor, '{Spec:1-Ax:X}Mtr')
    tz = Cpt(EpicsMotor, '{Spec:1-Ax:Z}Mtr')
    sy = Cpt(EpicsMotor, '{Env:1-Ax:Y}Mtr')
    sx = Cpt(EpicsMotor, '{Env:1-Ax:X}Mtr')
    sz = Cpt(EpicsMotor, '{Env:1-Ax:Z}Mtr')



analyzer = Analyzer('XF:10IDD', name='analyzer')
analyzer_xtals = AnalyzerDXtals('XF:10IDD-ES{Ez4:', name='analyzer_xtals')
analyzer_slits = AnalyzerSlits('XF:10IDD-ES{Ez4:', name='analyzer_slits')
mcm_slits = MCMSlits('XF:10IDD-OP{Ez4:', name='mcm_slits')
sample_stage = SampleStage('XF:10IDD-OP', name='sample_stage')

