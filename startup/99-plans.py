import bluesky.plans as bp
import bluesky.plan_stubs as bps
import bluesky.preprocessors as bpp
import bluesky.callbacks.fitting
import numpy as np
import lmfit
from bluesky.callbacks import LiveFit
from bluesky.suspenders import SuspendFloor
from ophyd import EpicsSignal
from tabulate import tabulate

tm1sum = EpicsSignal('XF:10ID-BI:TM176:SumAll:MeanValue_RBV')
susp = SuspendFloor(tm1sum, 1.e-5, resume_thresh = 1.e-5, sleep = 1*60)

def align_with_fit(dets, mtr, start, stop, gaps, md=None):
    # Performs relative scan of motor and retuns data staistics

    md = md or {}

    local_peaks = []
    for det in dets:
        for hint in det.hints['fields']:
            local_peaks.append(
                bluesky.callbacks.fitting.PeakStats(mtr.hints['fields'][0], hint)
            ) 
    # TODO use relative wrapper to avoid the reset behavior (or make it optional)
    plan = bpp.subs_wrapper(
        bp.rel_scan(dets, mtr, start, stop, gaps+1, md=md), 
        local_peaks
        )
    yield from plan
    return local_peaks

#def set_lambda_exposure(exposure):
#    det = lambda_det
#    yield from bps.mv(det.cam.acquire_time, exposure, det.cam.acquire_period, exposure)

def check_zero(dets=None, start=-20, stop=20, gaps=200, exp_time=1):
    # Performs relative scan of the HRM energy at tth = 0 and positions it to the peak center

    #ssxop = 0
    print('scanning zero')
    #yield from bps.mov(spec.tth, 0, spec.phi, 0, sample_stage.sx, -1)
    yield from bps.mv(spec.tth, 0)
    sample_pos = yield from bps.read(sample_stage)
    print(sample_pos)
    if dets is None:
        dets = [lambda_det]
        yield from set_lambda_exposure(exp_time)

    yield from bps.mv(whl, 7)
    for d in dets:
        # set the exposure times
        pass

    local_peaks = yield from align_with_fit(dets, hrmE, start, stop, gaps)
    
    cen = local_peaks[0].cen
    if cen is not None:
        target = 0.2 * (cen // .2)
        # move too far for backlash compensation
        yield from bps.mv(hrmE, target - 20)
        # apporach target from negative side 
        yield from bps.mv(hrmE, target)

def do_the_right_thing(i_time):
    yield from bps.mv(det1.integration_time, i_time)
    yield from count([det1])

def ct(exp_time):
    yield from bps.mv(sclr.preset_time, exp_time)
    yield from bp.count([sclr])


def double_ct(exp_time):
    yield from ct(exp_time)
    # yield from bps.mv(sample_stage.sx, 0)
    yield from ct(exp_time)

def Lipid_Qscan():
    # Test plan for the energy scan at several Q values
    tth001 = 16.8
    Qq = [1, 2, 3]
    c22 = sclr.channels.chan22
    yield from bps.mv(analyzer_slits.top, 1, analyzer_slits.bottom, -1, analyzer_slits.outboard, 1.5, analyzer_slits.inboard, -1.5)
    for kk in range(6):
        yield from bps.mv(anapd, 25)
        #yield from set_lambda_exposure(2)
        yield from check_zero(exp_time=2)
        yield from bps.mv(whl, 0)

        for q in Qq:
            th = qq2th(q)
            yield from bps.mv(spec.tth, th)
            yield from hrmE_dscan(-15, 15, 150, 30)

            yield from bps.mvr(sample_stage.sx, 0.03)
            yield from bps.mv(spec.tth, tth001)
            yield from set_lambda_exposure(5)
            loc_peaks = yield from align_with_fit([lambda_det], sample_stage.sy, -0.1, 0.1, 40)
            max_pos = local_peaks[0].max
            yield from bps.mvr(ample_stage.sy, -0.1)
            yield from bps.mv(ample_stage.sy, max_pos)
            loc_peaks = yield from align_with_fit([lambda_det], sample_stage.sz, -2, 2, 40)
            max_pos = local_peaks[0].max
            yield from bps.mv(sample_stage.sz, max_pos)
            loc_peaks = yield from align_with_fit([lambda_det], sample_stage.sy, -0.1, 0.1, 40)
            max_pos = local_peaks[0].max
            yield from bps.mvr(ample_stage.sy, -0.1)
            yield from bps.mv(ample_stage.sy, max_pos)

        yield from bps.mv(anapd, 3, spec.tth, 1)
        yield from bps.mv(sclr.channels.chan22.preset_time, 5)
        yield from bp.scan([c22], spec.tth, 1, 21, 101)

def Lipid_Qscan_wBC():
    # Lipid_Qscan with beam check
    yield from bpp.suspend_wrapper(Lipid_Qscan(), susp)


def GCarbon_Qscan():
    # Test plan for the energy resolution at Q=1.2 with the Glassy Carbon
    Qq = [1.2]
    yield from bps.mv(analyzer_slits.top, 1, analyzer_slits.bottom, -1, analyzer_slits.outboard, 1.5, analyzer_slits.inboard, -1.5)
    yield from bps.mv(anapd, 25, whl, 0)
    plt.clf()

    for kk in range(1):
        #yield from set_lambda_exposure(2)
     #   yield from check_zero(start=-10, stop=10, gaps=80,exp_time=2)
        for q in Qq:
            th = qq2th(q)
            yield from bps.mv(spec.tth, th)
            yield from hrmE_dscan(-10, 10, 100, 2)


def gaussian(x, A, sigma, x0):
    return A*np.exp(-(x - x0)**2/(2 * sigma**2))

#model = lmfit.Model(gaussian)
#init_guess = {'A': 800, 'sigma': 0.7, 'x0': 0}
#lf = LiveFit(model,'lambda_det_stats7_total', {'x': 'hrmE'}, init_guess)

def calc_lmfit(uid=-1, x="hrmE", channel=7):
    # Calculates fitting parameters for Gaussian function for energy scan with UID and Lambda channel
    hdr = db[uid]
    table = hdr.table()
    y = f'lambda_det_stats{channel}_total'
    lf = LiveFit(model, y, {'x': x}, {'A': table[y].max(), 'sigma': 0.7, 'x0': table[x][table[y].argmax()+1]})
    for name, doc in hdr.documents():
        lf(name, doc)
    gauss = gaussian(table[x], **lf.result.values)
    plt.plot(table[x], table[y], label=f"raw, channel={channel}", marker = 'o', linestyle = 'none')
    plt.plot(table[x], gauss.values, label=f"gaussian fit {channel}")
    plt.legend()
    return lf.result.values


def DxtalTempCalc(uid=-1):
    # Calculates temperature correction for the D crystals
    E0 = 9131.7     # energy (eV)
    TH = 88.5       # Dxtal asymmetry angle (deg)
    C1 = 3.725e-6   # constant (1/K)
    C2 = 5.88e-3    # constant (1/K)
    C3 = 5.548e-10  # constant (1/K2)
    T1 = 124.0      # temperature (K)
    T0 = 300.15     # crystal average temperature (K)

    Dtemp1 = EpicsSignal("XF:10ID-CT{FbPid:01}PID.VAL", name="Dtemp1")
    Dtemp2 = EpicsSignal("XF:10ID-CT{FbPid:02}PID.VAL", name="Dtemp2")
    Dtemp3 = EpicsSignal("XF:10ID-CT{FbPid:03}PID.VAL", name="Dtemp3")
    Dtemp4 = EpicsSignal("XF:10ID-CT{FbPid:04}PID.VAL", name="Dtemp4")
    Dtemp5 = EpicsSignal("XF:10ID-CT{FbPid:05}PID.VAL", name="Dtemp5")
    Dtemp6 = EpicsSignal("XF:10ID-CT{FbPid:06}PID.VAL", name="Dtemp6")

    bet = C1*(1 - np.exp(-C2*(T0-T1))) + C3*T0
    dE = []
    plt.clf()
    for n in range(1,7):
        fit_par = calc_lmfit(uid, channel=n)
        if fit_par['A'] < 100:
            print('**********************************')
            print('         WARNING !')
            print('      Fitting Error')
            return
        
        dE.append(fit_par['x0'])

    dE = [x-dE[0] for x in dE]
    dTe = [1.e-3*x/E0/bet for x in dE]
    dTh = [1.e3*x*np.tan(np.radians(TH))/E0 for x in dE]
    
    DTe = [Dtemp1.read()['Dtemp1']['value']+dTe[0], 
           Dtemp2.read()['Dtemp2']['value']+dTe[1], 
           Dtemp3.read()['Dtemp3']['value']+dTe[2], 
           Dtemp4.read()['Dtemp4']['value']+dTe[3], 
           Dtemp5.read()['Dtemp5']['value']+dTe[4], 
           Dtemp6.read()['Dtemp6']['value']+dTe[5]]
    Dheader = [' ', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6']
    dE.insert(0,'dEnrg')
    dTe.insert(0,'dTemp')
    dTh.insert(0,'dThe')
    DTe.insert(0,'Dtemp')
    Ddata = [dE, dTh, dTe, DTe]
    print('---------------------------------------------------------------------')
    print(tabulate(Ddata, headers=Dheader, tablefmt='pipe', stralign='center', floatfmt='.4f'))
    print('---------------------------------------------------------------------\n')
    update_temp = input('Do you want to update the temperature (yes/no): ')
    if update_temp == 'yes':
#        d1 = Dtemp1.set(DTe[1])
#        d2 = Dtemp2.set(DTe[2])
#        d3 = Dtemp3.set(DTe[3])
#        d4 = Dtemp4.set(DTe[4])
#        d5 = Dtemp5.set(DTe[5])
#        d6 = Dtemp6.set(DTe[6])
#        wait(d1, d2, d3, d4, d5, d6)
        print('\n')
        print('The temperature is updated')
    else:
        print('\n')
        print('Update is canceled')
    return {'dEn':dE, 'dTem':dTe, 'dThe':dTh, 'DTem':DTe}