from ophyd.controls import EpicsSignal, EpicsMotor

# SR current
# CNT012 src  SRcur  SR:C03-BI{DCCT:1}I:Real-I
SRcur = EpicsSignal('SR:C03-BI{DCCT:1}I:Real-I', rw=False, name='SRcur')

# Undulator: this is an example

#epu1_gap = PVPositioner('XF:23ID-ID{EPU:1-Ax:Gap}Pos-SP',
#                        readback='XF:23ID-ID{EPU:1-Ax:Gap}Pos-I',
#                        stop='SR:C23-ID:G1A{EPU:1-Ax:Gap}-Mtr.STOP',
#                        stop_val=1,
#                        put_complete=True,
#                        name='epu1_gap')

############# Front End Slits (Primary Slits) ############
# MOT001 fst  FSlitTop  FE:C10A-OP{Slt:1-Ax:T}Mtr.VAL
fst = EpicsMotor('FE:C10A-OP{Slt:1-Ax:T}Mtr', name='fst')

# MOT002 fsb  FSlitBot  FE:C10A-OP{Slt:2-Ax:B}Mtr.VAL
fsb = EpicsMotor('FE:C10A-OP{Slt:2-Ax:B}Mtr', name='fsb')

# MOT003 fso  FSlitOut  FE:C10A-OP{Slt:1-Ax:O}Mtr.VAL
fso = EpicsMotor('FE:C10A-OP{Slt:1-Ax:O}Mtr', name='fso')

# MOT004 fsi  FSlitIn  FE:C10A-OP{Slt:2-Ax:I}Mtr.VAL
fsi = EpicsMotor('FE:C10A-OP{Slt:2-Ax:I}Mtr', name='fsi')

''' Example of how to configure the Frontend slit composite motions
fe_xc = PVPositioner('FE:C23A-OP{Slt:12-Ax:X}center',
                     readback='FE:C23A-OP{Slt:12-Ax:X}t2.D',
                     stop='FE:C23A-CT{MC:1}allstop',
                     stop_val=1, put_complete=True,
                     name='fe_xc')

fe_yc = PVPositioner('FE:C23A-OP{Slt:12-Ax:Y}center',
                     readback='FE:C23A-OP{Slt:12-Ax:Y}t2.D',
                     stop='FE:C23A-CT{MC:1}allstop',
                     stop_val=1,
                     put_complete=True,
                     name='fe_yc')

fe_xg = PVPositioner('FE:C23A-OP{Slt:12-Ax:X}size',
                     readback='FE:C23A-OP{Slt:12-Ax:X}t2.C',
                     stop='FE:C23A-CT{MC:1}allstop',
                     stop_val=1, put_complete=True,
                     name='fe_xg')

fe_yg = PVPositioner('FE:C23A-OP{Slt:12-Ax:Y}size',
                     readback='FE:C23A-OP{Slt:12-Ax:Y}t2.C',
                     stop='FE:C23A-CT{MC:1}allstop',
                     stop_val=1,
                     put_complete=True,
                     name='fe_yg')
'''
