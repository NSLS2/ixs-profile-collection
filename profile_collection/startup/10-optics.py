from ophyd.controls import EpicsMotor, PVPositioner

############# DCM ############# 
# MOT005 dcm_y  DCM_Y XF:10IDA-OP{Mono:DCM-Ax:Y}Mtr.VAL
dcm_y = EpicsMotor('XF:10IDA-OP{Mono:DCM-Ax:Y}Mtr', name='dcm_y') 
# MOT006 dcm_th  DCM_theta  XF:10IDA-OP{Mono:DCM-Ax:P}Mtr.VAL
dcm_th = EpicsMotor('XF:10IDA-OP{Mono:DCM-Ax:P}Mtr', name='dcm_th') 
# MOT007 dcm_z2  DCM_Z2  XF:10IDA-OP{Mono:DCM-Ax:Z2}Mtr.VAL
dcm_z2 = EpicsMotor('XF:10IDA-OP{Mono:DCM-Ax:Z2}Mtr', name='dcm_z2')
# MOT008 dcm_p1  DCM_P1  XF:10IDA-OP{Mono:DCM-Ax:P1}Mtr.VAL
dcm_p1 = EpicsMotor('XF:10IDA-OP{Mono:DCM-Ax:P1}Mtr', name='dcm_p1') 
# MOT009 dcm_r2  DCM_R2  XF:10IDA-OP{Mono:DCM-Ax:R2}Mtr.VAL
dcm_r2 = EpicsMotor('XF:10IDA-OP{Mono:DCM-Ax:R2}Mtr', name='dcm_r2') 
# MOT010 dcm_pf  DCM_PF  XF:10IDA-OP{Mono:DCM-Ax:PF}Mtr.VAL
dcm_pf = EpicsMotor('XF:10IDA-OP{Mono:DCM-Ax:PF}Mtr', name='dcm_pf')

############## Mirrors  ############# 
# VFM
# MOT027 vmux  VFM_USX  XF:10IDD-OP{VFM:1-Ax:UX}Mtr.VAL
vmux = EpicsMotor('XF:10IDD-OP{VFM:1-Ax:UX}Mtr', name='vmux') 
# MOT028 vmuy  VFM_USY  XF:10IDD-OP{VFM:1-Ax:UY}Mtr.VAL
vmuy = EpicsMotor('XF:10IDD-OP{VFM:1-Ax:UY}Mtr', name='vmuy') 
# MOT029 vmdx  VFM_DSX  XF:10IDD-OP{VFM:1-Ax:DX}Mtr.VAL
vmdx = EpicsMotor('XF:10IDD-OP{VFM:1-Ax:DX}Mtr', name='vmdx') 
# MOT030 vmdy  VFM_DSY  XF:10IDD-OP{VFM:1-Ax:DY}Mtr.VAL
vmdy = EpicsMotor('XF:10IDD-OP{VFM:1-Ax:DY}Mtr', name='vmdy')

# HFM
# MOT031 hmux  HFM_USX  XF:10IDD-OP{HFM:1-Ax:UX}Mtr.VAL
hmux = EpicsMotor('XF:10IDD-OP{HFM:1-Ax:UX}Mtr', name='hmux')
# MOT032 hmuy  HFM_USY  XF:10IDD-OP{HFM:1-Ax:UY}Mtr.VAL
hmuy = EpicsMotor('XF:10IDD-OP{HFM:1-Ax:UY}Mtr', name='hmuy')
# MOT033 hmdx  HFM_DSX  XF:10IDD-OP{HFM:1-Ax:DX}Mtr.VAL
hmdx = EpicsMotor('XF:10IDD-OP{HFM:1-Ax:DX}Mtr', name='hmdx')
# MOT034 hmdy  HFM_DSY  XF:10IDD-OP{HFM:1-Ax:DY}Mtr.VAL
hmdy = EpicsMotor('XF:10IDD-OP{HFM:1-Ax:DY}Mtr', name='hmdy')


############## Slits ############### 
# MOT017 s1t  Slit1Top  XF:10IDA-OP{Slt:1-Ax:T}Mtr.VAL
s1t = EpicsMotor('XF:10IDA-OP{Slt:1-Ax:T}Mtr', name='s1t')
# MOT018 s1b  Slit1Bot  XF:10IDA-OP{Slt:1-Ax:B}Mtr.VAL
s1b = EpicsMotor('XF:10IDA-OP{Slt:1-Ax:B}Mtr', name='s1b') 
# MOT019 s1o  Slit1Out  XF:10IDA-OP{Slt:1-Ax:O}Mtr.VAL
s1o = EpicsMotor('XF:10IDA-OP{Slt:1-Ax:O}Mtr', name='s1o') 
# MOT020 s1i  Slit1In  XF:10IDA-OP{Slt:1-Ax:I}Mtr.VAL
s1i = EpicsMotor('XF:10IDA-OP{Slt:1-Ax:I}Mtr', name='s1i')
# x width
s1xg = PVPositioner('XF:10IDA-OP{Slt:1-Ax:X}size',
                    readback='XF:10IDA-OP{Slt:1-Ax:X}t2.C',
                    done='XF:10IDA-OP{Slt:1-Ax:X}DMOV',
                    done_val=1,
                    name='s1xg')
s1xc = PVPositioner('XF:10IDA-OP{Slt:1-Ax:X}center',
                    readback='XF:10IDA-OP{Slt:1-Ax:X}t2.D',
                    done='XF:10IDA-OP{Slt:1-Ax:X}DMOV',
                    done_val=1,
                    name='s1xc')
s1yg = PVPositioner('XF:10IDA-OP{Slt:1-Ax:Y}size',
                    readback='XF:10IDA-OP{Slt:1-Ax:Y}t2.C',
                    done='XF:10IDA-OP{Slt:1-Ax:Y}DMOV',
                    done_val=1,
                    name='s1yg')
s1yc = PVPositioner('XF:10IDA-OP{Slt:1-Ax:Y}center',
                    readback='XF:10IDA-OP{Slt:1-Ax:Y}t2.D',
                    done='XF:10IDA-OP{Slt:1-Ax:Y}DMOV',
                    done_val=1,
                    name='s1xy')

##
# SSA
# MOT025 ssat  SSA_T  XF:10IDB-OP{SSA:1-Ax:T}Mtr.VAL
ssat = EpicsMotor('XF:10IDB-OP{SSA:1-Ax:T}Mtr', name='ssat')
# MOT026 ssab  SSA_B  XF:10IDB-OP{SSA:1-Ax:B}Mtr.VAL
ssab = EpicsMotor('XF:10IDB-OP{SSA:1-Ax:B}Mtr', name='ssab')

# MOT021 s2t  Slit2Top  XF:10IDC-OP{Slt:4-Ax:T}Mtr.VAL
s2t = EpicsMotor('XF:10IDC-OP{Slt:4-Ax:T}Mtr', name='s2t') 
# MOT022 s2b  Slit2Bot  XF:10IDC-OP{Slt:4-Ax:B}Mtr.VAL
s2b = EpicsMotor('XF:10IDC-OP{Slt:4-Ax:B}Mtr', name='s2b') 
# MOT023 s2o  Slit2Out  XF:10IDC-OP{Slt:4-Ax:O}Mtr.VAL
s2o = EpicsMotor('XF:10IDC-OP{Slt:4-Ax:O}Mtr', name='s2o') 
# MOT024 s2i  Slit2In  XF:10IDC-OP{Slt:4-Ax:I}Mtr.VAL
s2i = EpicsMotor('XF:10IDC-OP{Slt:4-Ax:I}Mtr', name='s2i')

############## Diagnostic Manipulators ############## 
# MOT011 bpm1_y  BPM1_Y  XF:10IDA-BI{BPM:1-Ax:YFoil}Mtr.VAL
bpm1_y = EpicsMotor('XF:10IDA-BI{BPM:1-Ax:YFoil}Mtr', name='bpm1_y')
# MOT012 bpm1_dx  BPM1_DX  XF:10IDA-OP{BPM:1-Ax:X}Mtr.VAL
bpm1_dx = EpicsMotor('XF:10IDA-OP{BPM:1-Ax:X}Mtr', name='bpm1_dx') 
# MOT013 bpm1_dy  BPM1_DY  XF:10IDA-OP{BPM:1-Ax:Y}Mtr.VAL
bpm1_dy = EpicsMotor('XF:10IDA-OP{BPM:1-Ax:Y}Mtr', name='bpm1_dy') 
# MOT014 bpm2_y  BPM2_Y  XF:10IDC-BI{BPM:2-Ax:Y}Mtr.VAL
bpm2_y = EpicsMotor('XF:10IDC-BI{BPM:2-Ax:Y}Mtr', name='bpm2_y') 
# MOT015 bpm2_dx  BPM2_DX  XF:10IDC-OP{BPM:2-Ax:X}Mtr.VAL
bpm2_dx = EpicsMotor('XF:10IDC-OP{BPM:2-Ax:X}Mtr', name='bpm2_dx') 
# MOT016 bpm2_dy  BPM2_DY  XF:10IDC-OP{BPM:2-Ax:Y}Mtr.VAL
bpm2_dy = EpicsMotor('XF:10IDC-OP{BPM:2-Ax:Y}Mtr', name='bpm2_dy') 
