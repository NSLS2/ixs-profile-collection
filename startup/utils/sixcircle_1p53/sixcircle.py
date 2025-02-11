# sixcircle: code for six-circle diffractometer
# Required additional documents: 1- scbasic.py, 2- ini.conf

# Developed with Python 3.7.3, numpy 1.16.4, scipy 1.3.0
# Materials Dynamics Laboratory, RIKEN SPring-8 Center
# Ver 1.51, December 2020

# Authors: Wenyang Zhao & Alfred Q. R. Baron
# Contact: baron@spring8.or.jp

# If this code is used independently of work at BL43, please reference the writeup ("Open-source Python software for six-circle diffraction with an inelastic x-ray scattering (IXS) spectrometer" by Wenyang ZHAO and Alfred Q.R. BARON, unpublished, available at https://beamline.harima.riken.jp/bl43lxu)

# In typical use related to experimental work at BL43LXU, it is enough to reference a BL publication.  However, if this code is used extensively, then please include a separate reference as above.

import numpy as np

import os.path
from os import path

print ('')
print ('RSC Materials Dynamics Laboratory Six Circle Code by Wenyang Zhao and Alfred Q.R. Baron')
print ('')
print ('Suggest:  import sixcircle; from sixcircle import *')
print ('For IXS users at SPring-8:  from sixcircle_rqd import *')
print ('')

# Module of six-circle calculation
from . import scbasic

# allow some output to screen to be skipped
runquiet = False
def runquiet_on():
    global runquiet
    runquiet = True
def runquiet_off():
    global queitrunning
    runquiet = False

def wdesc(fname,dstr):
    print(' - {:40}'.format(fname),end='')
    print(' {}'.format(dstr))

print('Function Definitions from sixcicle.py:')
wdesc("ini()","SixCircle initialization")
# Initialization
def ini():
    print('\nInitilizing sixcircle setup.')
    # Valid mnemonics
    global mnemonics
    mnemonics = ['tth','th','chi','phi','mu','gam']
    # Positions: all zero at startup
    global TTH, TH, CHI, PHI, MU, GAM
    TTH=0.0; TH=0.0; CHI=0.0; PHI=0.0; MU=0.0; GAM=0.0
    # SA: scattering angle; ABSQ: length of scattering vector
    global SA, ABSQ
    SA = 0.0; ABSQ = 0.0
    # Omega: difference from TTH/2 to TH
    global OMEGA
    OMEGA = TH - TTH/2
    # Frozen of six-circle calculation. Default '456' at startup (freeze mu, gam, omega)
    global g_frozen
    g_frozen = '456'
    # Frozen positions: all zero at startup
    global F_TTH, F_TH, F_CHI, F_PHI, F_MU, F_GAM, F_OMEGA, F_AZIMUTH, F_ALPHA, F_BETA
    F_TTH=0.0; F_TH=0.0; F_CHI=0.0; F_PHI=0.0; F_MU=0.0; F_GAM=0.0; F_OMEGA=0.0; F_AZIMUTH=0.0; F_ALPHA=0.0; F_BETA=0.0
    # Output precision, default PRE=3, (PRE+2) for wavelength
    global PRE
    PRE = 3
    # Flag of using wh() after mv() or br(), default False
    global FLAG_WH
    FLAG_WH = True
    # Load default configuration file at startup
    
    if path.isfile('/nsls2/data3/ixs/shared/config/bluesky/profile_collection/startup/utils/sixcircle_1p53/sixcircle_last_UB') :
        load('/nsls2/data3/ixs/shared/config/bluesky/profile_collection/startup/utils/sixcircle_1p53/sixcircle_last_UB')
    else :
        load('/nsls2/data3/ixs/shared/config/bluesky/profile_collection/startup/utils/sixcircle_1p53/ini.conf')

# Load a configuration file
wdesc ("load('filepath'), save('filepath')","Loads,Saves Wavelength and orientation matrix")
def load(filepath):
    # Sample description
    global g_sample
    # Azimuth reference vector H, K, L
    global g_haz, g_kaz, g_laz
    # Lattice parameters: a, b, c, alpha, beta, gam
    global g_aa, g_bb, g_cc, g_al, g_be, g_ga
    # Primary reflection: H, K, L
    global g_h0, g_k0, g_l0
    # Primary reflection: positions of tth, th, chi, phi, mu, gam
    global g_u00, g_u01, g_u02, g_u03, g_u04, g_u05
    # Secondary reflections: H, K, L
    global g_h1, g_k1, g_l1
    # Secondary reflection: positions of tth, th, chi, phi, mu, gam 
    global g_u10, g_u11, g_u12, g_u13, g_u14, g_u15
    # Wavelength in finding reference reflections
    global g_lambda0, g_lambda1
    # Wavelength in current calculation
    global LAMBDA
    # Limit of positions
    global L_TTH, U_TTH, L_TH, U_TH, L_CHI, U_CHI, L_PHI, U_PHI, L_MU, U_MU, L_GAM, U_GAM, L_ALPHA, U_ALPHA, L_BETA, U_BETA
    try:
        with open(filepath, 'r') as f:
            content = f.readlines()
        if not runquiet : print ('Reading configuration from {0}'.format(filepath))
    except:
        print ('\nError reading configuration file {0}'.format(filepath))
        return
    dic = {}
    for line in content:
            if line.startswith('GLOBAL') == True:
                gvar = line.split()[1]
                gvarvalue = line.split()[2]
                dic.update({gvar:gvarvalue})
    g_sample = dic['g_sample']
    g_haz = float(dic['g_haz'])
    g_kaz = float(dic['g_kaz'])
    g_laz = float(dic['g_laz'])
    g_aa = float(dic['g_aa'])
    g_bb = float(dic['g_bb'])
    g_cc = float(dic['g_cc'])
    g_al = float(dic['g_al'])
    g_be = float(dic['g_be'])
    g_ga = float(dic['g_ga'])
    g_lambda0 = float(dic['g_lambda0'])
    g_h0 = float(dic['g_h0'])
    g_k0 = float(dic['g_k0'])
    g_l0 = float(dic['g_l0'])
    g_u00 = float(dic['g_u00'])
    g_u01 = float(dic['g_u01'])
    g_u02 = float(dic['g_u02'])
    g_u03 = float(dic['g_u03'])
    g_u04 = float(dic['g_u04'])
    g_u05 = float(dic['g_u05'])
    g_lambda1 = float(dic['g_lambda1'])
    g_h1 = float(dic['g_h1'])
    g_k1 = float(dic['g_k1'])
    g_l1 = float(dic['g_l1'])
    g_u10 = float(dic['g_u10'])
    g_u11 = float(dic['g_u11'])
    g_u12 = float(dic['g_u12'])
    g_u13 = float(dic['g_u13'])
    g_u14 = float(dic['g_u14'])
    g_u15 = float(dic['g_u15'])
    L_TTH = float(dic['L_TTH'])
    U_TTH = float(dic['U_TTH'])
    L_TH = float(dic['L_TH'])
    U_TH = float(dic['U_TH'])
    L_CHI = float(dic['L_CHI'])
    U_CHI = float(dic['U_CHI'])
    L_PHI = float(dic['L_PHI'])
    U_PHI = float(dic['U_PHI'])
    L_MU = float(dic['L_MU'])
    U_MU = float(dic['U_MU'])
    L_GAM = float(dic['L_GAM'])
    U_GAM = float(dic['U_GAM'])
    L_ALPHA = float(dic['L_ALPHA'])
    U_ALPHA = float(dic['U_ALPHA'])
    L_BETA = float(dic['L_BETA'])
    U_BETA = float(dic['U_BETA'])
    # Set wavelength in current calculation
    LAMBDA = g_lambda0
    # Calculate matrix U and B
    UB()

# Save a configuration file
#wdesc ("save('filepath')","Saves crystal parameters, Wavelength, and orientation matrix")
def save(filepath):
    try:
        with open(filepath, 'w') as f:
            f.write('# Configuration file of sixcircle.\n')
            f.write('\n')
            f.write('# Sample description\n')
            f.write('GLOBAL g_sample {0}\n'.format(g_sample))
            f.write('\n')
            f.write('# Azimuthal reference H K L\n')
            dic = dict(g_haz=g_haz, g_kaz=g_kaz, g_laz=g_laz)
            for key in dic.keys():
                f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],PRE))
            f.write('\n')
            f.write('# Lattice parameters\n')
            dic = dict(g_aa=g_aa, g_bb=g_bb, g_cc=g_cc, g_al=g_al, g_be=g_be, g_ga=g_ga)
            for key in dic.keys():
                f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],PRE))
            f.write('\n')
            f.write('# Primary-reflection wavelength\n')
            f.write('GLOBAL g_lambda0 {0:.{1}f}\n'.format(g_lambda0,PRE+2))
            f.write('\n')
            f.write('# Primary-reflection HKL coordinates\n')
            dic = dict(g_h0=g_h0, g_k0=g_k0, g_l0=g_l0)
            for key in dic.keys():
                f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],PRE))
            f.write('\n')
            f.write('# Primary-reflection angles: tth, th, chi, phi, mu, gam\n')
            dic = dict(g_u00=g_u00, g_u01=g_u01, g_u02=g_u02, g_u03=g_u03, g_u04=g_u04, g_u05=g_u05)
            for key in dic.keys():
                f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],PRE))
            f.write('\n')
            f.write('# Secondary-reflection wavelength\n')
            f.write('GLOBAL g_lambda1 {0:.{1}f}\n'.format(g_lambda1,PRE+2))
            f.write('\n')
            f.write('# Secondary-reflection HKL coordinates\n')
            dic = dict(g_h1=g_h1, g_k1=g_k1, g_l1=g_l1)
            for key in dic.keys():
                f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],PRE))
            f.write('\n')
            f.write('# Secondary-reflection angles: tth, th, chi, phi, mu, gam\n')
            dic = dict(g_u10=g_u10, g_u11=g_u11, g_u12=g_u12, g_u13=g_u13, g_u14=g_u14, g_u15=g_u15)
            for key in dic.keys():
                f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],PRE))
            f.write('\n')
            f.write('# Limit of positions\n')
            dic1 = dict(L_TTH=L_TTH, U_TTH=U_TTH, L_TH=L_TH, U_TH=U_TH, L_CHI=L_CHI, U_CHI=U_CHI, L_PHI=L_PHI, U_PHI=U_PHI)
            dic2 = dict(L_MU=L_MU, U_MU=U_MU, L_GAM=L_GAM, U_GAM=U_GAM, L_ALPHA=L_ALPHA, U_ALPHA=U_ALPHA, L_BETA=L_BETA, U_BETA=U_BETA)
            for key in dic1.keys():
                f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic1[key],PRE))
            for key in dic2.keys():
                f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic2[key],PRE))
        if not runquiet : print ('\nWrote configuration to {0}'.format(filepath))
    except:
        print ('\nError in writing configuration file {0}'.format(filepath))

wdesc("UB()","Updates UB matrix")
def UB():
    # Lattice parameters in reciprocal space: a, b, c, alpha, beta, gam
    global g_aa_s, g_bb_s, g_cc_s, g_al_s, g_be_s, g_ga_s
    # Matrix B and matrix U
    global M_B, M_U
    # wavelength in finding reference reflections
    tupleB = scbasic.B_matrix(g_aa, g_bb, g_cc, g_al, g_be, g_ga)
    if (tupleB[0]) == False:
        print ('\nError! Invalid lattice parameters!')
        return
    flag, M_B, g_aa_s, g_bb_s, g_cc_s, g_al_s, g_be_s, g_ga_s = tupleB
    # Primary and secondary reflections: unit vectors calcualted from positions
    u0phi = scbasic.uphi_vector(g_u00/2, g_u01-g_u00/2, g_u02, g_u03, g_u04, g_u05)
    u1phi = scbasic.uphi_vector(g_u10/2, g_u11-g_u10/2, g_u12, g_u13, g_u14, g_u15)
    # Primary and secondary reflections: unit vectors cacluated from H, K, L
    h0b = np.array([[g_h0],[g_k0],[g_l0]],dtype=float)
    h1b = np.array([[g_h1],[g_k1],[g_l1]],dtype=float)
    u0c = scbasic.uc_vector(M_B, h0b)
    u1c = scbasic.uc_vector(M_B, h1b)
    # Calculate orientation matrix U
    errcode, M_U = scbasic.U_matrix(u0phi, u1phi, u0c, u1c)
    if errcode == 0:
#        print ('\nUB recalculated for {} '.format(g_sample),end='')
#        print (' with ({0:.{3}f}, {1:.{3}f}, {2:.{3}f} '.format(g_aa,g_bb,g_cc,PRE),end='')
#        print (' {1:.{3}f}, {2:.{3}f}, {2:.{3}f})'.format(g_al,g_be,g_ga,2))
#        print ('     or0 = ({0:.{3}f}, {1:.{3}f}, {2:.{3}f})    and '.format(g_h0,g_k0,g_l0,PRE),end='')
#        print ('     or1 = ({0:.{3}f}, {1:.{3}f}, {2:.{3}f}) '.format(g_h1,g_k1,g_l1,PRE))
        # Get (H, K, L, ALPHA, BETA, AZIMUTH) using new orientation matrix U
        wh_refresh()
        or_check()
        runquiet_on()
        save('sixcircle_last_UB')
        runquiet_off()
    elif errcode == 1:
        print ('\nCannot find orientation matrix:  Reflections are parallel.\n')
    elif errcode == 2:
        print ('\nCannot find orientation matrix:  Reflections (by angles) are parallel.\n')

# Show parameters of orientation calculation
wdesc ("pa()","Display crystal structure and reference vectors")
def pa():
    print ('')
    print ('Primary Reflection (or0, at lambda {0:.{1}f}):'.format(g_lambda0,PRE+2))
    print ('           tth={0:.{6}f}, th={1:.{6}f}, chi={2:.{6}f}, phi={3:.{6}f}, mu={4:.{6}f}, gam={5:.{6}f}'.format(g_u00,g_u01,g_u02,g_u03,g_u04,g_u05,PRE))
    print ('           H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(g_h0,g_k0,g_l0,PRE))
    print ('')
    print ('Secondary Reflection (or1, at lambda {0:.{1}f}):'.format(g_lambda1,PRE+2))
    print ('           tth={0:.{6}f}, th={1:.{6}f}, chi={2:.{6}f}, phi={3:.{6}f}, mu={4:.{6}f}, gam={5:.{6}f}'.format(g_u10,g_u11,g_u12,g_u13,g_u14,g_u15,PRE))
    print ('           H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(g_h1,g_k1,g_l1,PRE))
    print ('')
    print ('Lattice Constants (lengths / angles):')
    print ('           real space = {0:.{6}f} {1:.{6}f} {2:.{6}f} / {3:.{6}f} {4:.{6}f} {5:.{6}f}'.format(g_aa,g_bb,g_cc,g_al,g_be,g_ga,PRE))
    print ('           reciprocal space = {0:.{6}f} {1:.{6}f} {2:.{6}f} / {3:.{6}f} {4:.{6}f} {5:.{6}f}'.format(g_aa_s,g_bb_s,g_cc_s,g_al_s,g_be_s,g_ga_s,PRE))
    print ('')
    print ('Azimuthal Reference:')
    print ('           H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(g_haz, g_kaz, g_laz,PRE))
    print ('')
    print ('           LAMBDA = {0:.{1}f}'.format(LAMBDA,PRE+2))
    print ('')

# Set wavelength
wdesc('setlambda(wavelength)','Sets LAMBDA in angstroms')
def setlambda(*args):
    global LAMBDA
    if len(args) == 0:
        LAMBDA_input =  input('\nWavelength / A ({0:.{1}f})? '.format(LAMBDA,PRE+2))
        if len(LAMBDA_input) != 0:
            LAMBDA = float(LAMBDA_input)
    elif len(args) == 1:
        if (type(args[0]) in [int,float]) == False:
            print ('\nInvalid argument: {0}\n'.format(args[0]))
            return
        LAMBDA = float(args[0])
    else:
        print ('\nUsage:  setlambda()  or  setlambda(LAMBDA)\n')
        return
    print("  -> LAMBDA now %.8f A" %(LAMBDA))
    wh_refresh()

# Set lattice parameters
# Usage: setlat() or setlat(a,b,c,alpha,beta,gam)
wdesc('setlat(a,b,c,alpha,beta,gamma)','Set crystal parameters (A & deg)')
def setlat(*args):
    global g_sample, g_aa, g_bb, g_cc, g_al, g_be, g_gl, g_be, g_ga
    if len(args) == 0:
        print ('\nEnter real space lattice parameters:')
        g_aa_input = input(' Lattice a ({0:.{1}f})? '.format(g_aa,PRE))
        if len(g_aa_input) != 0:
            g_aa = float(g_aa_input)
        g_bb_input = input(' Lattice b ({0:.{1}f})? '.format(g_bb,PRE))
        if len(g_bb_input) != 0:
            g_bb = float(g_bb_input)
        g_cc_input = input(' Lattice c ({0:.{1}f})? '.format(g_cc,PRE))
        if len(g_cc_input) != 0:
            g_cc = float(g_cc_input)
        g_al_input = input(' Lattice alpha ({0:.{1}f})? '.format(g_al,PRE))
        if len(g_al_input) != 0:
            g_al = float(g_al_input)
        g_be_input = input(' Lattice beta ({0:.{1}f})? '.format(g_be,PRE))
        if len(g_be_input) != 0:
            g_be = float(g_be_input)
        g_ga_input = input(' Lattice gam ({0:.{1}f})? '.format(g_ga,PRE))
        if len(g_ga_input) != 0:
            g_ga = float(g_ga_input)
    elif len(args) == 6:
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        g_aa, g_bb, g_cc, g_al, g_be, g_ga = args
    else:
        print ('\nUsage:  setlat()  or  setlat(a,b,c,alpha,beta,gam)\n')
        return
    g_sample_input = input('\nSample description: ({0})? '.format(g_sample))
    if len(g_sample_input) != 0:
        g_sample = g_sample_input
    print ('    -> Sample name set to {0}'.format(g_sample))
    # Calculate matrix U and B using new lattice parameters
    UB()
    # Get (H, K, L, ALPHA, BETA, AZIMUTH) using new lattice parameters
    wh_refresh()

# Set frozen of six-circle calculation
# Usage: setfrozen() or e.g., setfrozen(456)  or  setmode('045')
wdesc("setfrozen() , setfrozen(456)","Choose which angles to freeze")
def setfrozen(*args):
    global g_frozen
    dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
    dic_pos = {0:F_TTH, 1:F_TH, 2:F_CHI, 3:F_PHI, 4:F_MU, 5:F_GAM, 6:F_OMEGA, 7:F_AZIMUTH, 8:F_ALPHA, 9:F_BETA}
    if len(args) == 0:
        print ('')
        print ('Current frozen: {0}'.format(g_frozen))
        g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(g_frozen)]
        print ('Current frozen angles:')
        print (' {0:>20}{1:>20}{2:>20}'.format(dic_ang[g_frozen_1],dic_ang[g_frozen_2],dic_ang[g_frozen_3]))
        print (' {0:>20.{3}f}{1:>20.{3}f}{2:>20.{3}f}'.format(dic_pos[g_frozen_1],dic_pos[g_frozen_2],dic_pos[g_frozen_3],PRE))
        print ('')
        print ('tth(0)  th(1)  chi(2)  phi(3)  mu(4)  gam(5)  omega(6)  azimuth(7)  alpha(8)  beta(9)')
        loop_flag = True
        while loop_flag == True:
            print ('')
            g_frozen_input = input('Select three frozen angles (a three-digit integer, e.g. 456): ')
            # If no input, use current setfrozen
            if len(g_frozen_input) == 0:
                g_frozen_input = str(g_frozen)
                loop_flag = False
            # Check input is a three-digit integer
            if g_frozen_input.isdigit()==False or len(g_frozen_input)!=3:
                print ('Invalid argument: not a three-digit integer.')
                continue
            # Check input has repeated number
            rep_flag = False
            for i in range(0,10):
                if g_frozen_input.count(str(i)) >= 2:
                    rep_flag = True
                    continue
            if rep_flag == True:
                print ('Invalid argument: due to repetition.')
                continue
            # Check at most one frozen angle in {tth, th, omega}
            if g_frozen_input.count('0') + g_frozen_input.count('1') + g_frozen_input.count('6') > 1:
                print ('Invalid frozen: at most one frozen angle in {tth, th, omega}')
                continue
            # Check at most two frozen angles in {tth, mu, gam}
            if g_frozen_input.count('0') + g_frozen_input.count('4') + g_frozen_input.count('5') > 2:
                print ('Invalid frozen: at most two frozen angles in {tth, mu, gam}')
                continue
            # Check at most one frozen angle in {azimuth, alpha, beta}
            if g_frozen_input.count('7') + g_frozen_input.count('8') + g_frozen_input.count('9') > 1:
                print ('Invalid frozen: at most one frozen angle in {azimuth, alpha, beta}')
                continue
            # Valid input
            loop_flag = False
    if len(args) == 1:
        g_frozen_input = str(args[0])
        if g_frozen_input.isdigit()==False or len(g_frozen_input)!=3:
            print ('Invalid argument: not a three-digit integer.')
            return
        for i in range(0,10):
            if g_frozen_input.count(str(i)) >= 2:
                print ('Invalid argument: due to repetition.')
                return
        if g_frozen_input.count('0') + g_frozen_input.count('1') + g_frozen_input.count('6') > 1:
            print ('Invalid frozen: at most one frozen angle in {tth, th, omega}')
            return
        if g_frozen_input.count('0') + g_frozen_input.count('4') + g_frozen_input.count('5') > 2:
            print ('Invalid frozen: at most two fixed angles in {tth, mu, gam}')
            return
        if g_frozen_input.count('7') + g_frozen_input.count('8') + g_frozen_input.count('9') > 1:
            print ('Invalid frozen: at most one fixed angle in {azimuth, alpha, beta}')
            return
    if len(args) != 0 and len(args) != 1:
        print ("\nUsage:  setfrozen()  or  e.g., setfrozen('456')\n")
        return
    # New frozen
    g_frozen = ''.join(sorted(list(g_frozen_input)))
    g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(g_frozen)]
    print ('')
    print ('Current frozen: {0}'.format(g_frozen))
    print ('Current frozen angles:')
    print (' {0:>20}{1:>20}{2:>20}'.format(dic_ang[g_frozen_1],dic_ang[g_frozen_2],dic_ang[g_frozen_3]))
    print (' {0:>20.{3}f}{1:>20.{3}f}{2:>20.{3}f}'.format(dic_pos[g_frozen_1],dic_pos[g_frozen_2],dic_pos[g_frozen_3],PRE))
    print ('Use freeze() command to change frozen values.')
    print ('')

# Set positions of three frozen angles
# Usage: freeze() or freeze(position1, position2, position3)
wdesc ('freeze(), freeze(a1,a2,a3)','Choose values for frozen angles (degrees)')
def freeze(*args):
    dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
    dic_pos = {0:F_TTH, 1:F_TH, 2:F_CHI, 3:F_PHI, 4:F_MU, 5:F_GAM, 6:F_OMEGA, 7:F_AZIMUTH, 8:F_ALPHA, 9:F_BETA}
    g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(g_frozen)]
    def set_freeze_positions(g_frozen_i,position_i):
        if g_frozen_i == 0:
            global F_TTH
            F_TTH = position_i
        elif g_frozen_i == 1:
            global F_TH
            F_TH = position_i
        elif g_frozen_i == 2:
            global F_CHI
            F_CHI = position_i
        elif g_frozen_i == 3:
            global F_PHI
            F_PHI = position_i
        elif g_frozen_i == 4:
            global F_MU
            F_MU = position_i
        elif g_frozen_i == 5:
            global F_GAM
            F_GAM = position_i
        elif g_frozen_i == 6:
            global F_OMEGA
            F_OMEGA = position_i
        elif g_frozen_i == 7:
            global F_AZIMUTH
            F_AZIMUTH = position_i
        elif g_frozen_i == 8:
            global F_ALPHA
            F_ALPHA = position_i
        elif g_frozen_i == 9:
            global F_BETA
            F_BETA = position_i
    if len(args) == 0:
        print('')
        position_1_input = input(' Freeze {0} ({1:.{2}f})? '.format(dic_ang[g_frozen_1], dic_pos[g_frozen_1], PRE))
        if len(position_1_input) != 0:
            position_1 = float(position_1_input)
            set_freeze_positions(g_frozen_1,position_1)
        position_2_input = input(' Freeze {0} ({1:.{2}f})? '.format(dic_ang[g_frozen_2], dic_pos[g_frozen_2], PRE))
        if len(position_2_input) != 0:
            position_2 = float(position_2_input)
            set_freeze_positions(g_frozen_2,position_2)
        position_3_input = input(' Freeze {0} ({1:.{2}f})? '.format(dic_ang[g_frozen_3], dic_pos[g_frozen_3], PRE))
        if len(position_3_input) != 0:
            position_3 = float(position_3_input)
            set_freeze_positions(g_frozen_3,position_3)
    elif len(args) == 3:
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        position_1,position_2,position_3 = args
        set_freeze_positions(g_frozen_1,position_1)
        set_freeze_positions(g_frozen_2,position_2)
        set_freeze_positions(g_frozen_3,position_3)
    else:
        print ('\nUsage:  freeze()  or  freeze(position1, position2, position3)\n')
        return
    # Update positions of frozen angles
    dic_pos = {0:F_TTH, 1:F_TH, 2:F_CHI, 3:F_PHI, 4:F_MU, 5:F_GAM, 6:F_OMEGA, 7:F_AZIMUTH, 8:F_ALPHA, 9:F_BETA}
    print ('')
    print ('Positions of frozen angles:')
    print ('{0:>10}{1:>10}{2:>10}'.format(dic_ang[g_frozen_1],dic_ang[g_frozen_2],dic_ang[g_frozen_3]))
    print ('{0:>10.{3}f}{1:>10.{3}f}{2:>10.{3}f}'.format(dic_pos[g_frozen_1],dic_pos[g_frozen_2],dic_pos[g_frozen_3],PRE))

# Set azimuth reference vector
# Usage: setaz() or setaz(H,K,L)
wdesc('setaz(H,K,L)','Sets surface normal/azimuthal reference')
def setaz(*args):
    global g_haz, g_kaz, g_laz
    if len(args) == 0:
        print ('\nEnter azimuthal reference H K L:')
        g_haz_input = input(' Azimuthal H ({0})? '.format(g_haz))
        if len(g_haz_input) != 0:
            g_haz = float(g_haz_input)
        g_kaz_input = input(' Azimuthal K ({0})? '.format(g_kaz))
        if len(g_kaz_input) != 0:
            g_kaz = float(g_kaz_input)        
        g_laz_input = input(' Azimuthal L ({0})? '.format(g_laz))
        if len(g_laz_input) != 0:
            g_laz = float(g_laz_input)
    elif len(args) == 3:
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        g_haz, g_kaz, g_laz = args
    else:
        print ('\nUsage:  setaz()  or  setaz(H,K,L)\n')
        return
    # Get (..., ALPHA, BETA, AZIMUTH) using new azimuth reference vector
    wh_refresh()

# Set H, K, L of primary reflection
# Positions of primary reflection are current positions
# Usage: or0() or or0(H,K,L)
wdesc('or0(H,K,L) , or1(H,K,L)','Set primary,secondary (or0,or1) at present angle values')
def or0(*args):
    global g_u00, g_u01, g_u02, g_u03, g_u04, g_u05
    global g_h0, g_k0, g_l0
    global g_lambda0
    g_lambda0 = LAMBDA
    g_u00, g_u01, g_u02, g_u03, g_u04, g_u05 = (TTH, TH, CHI, PHI, MU, GAM)
    if len(args) == 0:
        print ('\nEnter primary-reflection HKL coordinates:')
        g_h0_input = input(' H ({0:.{1}f})? '.format(g_h0,PRE))
        g_k0_input = input(' K ({0:.{1}f})? '.format(g_k0,PRE))
        g_l0_input = input(' L ({0:.{1}f})? '.format(g_l0,PRE))
        if len(g_h0_input) != 0:
            g_h0 = float(g_h0_input)
        if len(g_k0_input) != 0:
            g_k0 = float(g_k0_input)
        if len(g_l0_input) != 0:
            g_l0 = float(g_l0_input)
    elif len(args) == 3:
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        g_h0, g_k0, g_l0 = args
    else:
        print ('\nUsage:  or0()  or  or0(H,K,L)\n')
        return
    UB()
    wh_refresh()

# Set H, K, L and positions of primary reflection
wdesc('setor0() , setor1()','Set primary,secondary (or0,or1) at entered angles')
def setor0():
    global g_u00, g_u01, g_u02, g_u03, g_u04, g_u05
    global g_h0, g_k0, g_l0
    global g_lambda0
    g_lambda0 = LAMBDA
    print ('\nEnter primary-reflection angles:')
    g_u00_input = input(' Two Theta ({0:.{1}f})? '.format(g_u00,PRE))
    if len(g_u00_input) != 0:
        g_u00 = float(g_u00_input)
    g_u01_input = input(' Theta ({0:.{1}f})? '.format(g_u01,PRE))
    if len(g_u01_input) != 0:
        g_u01 = float(g_u01_input)
    g_u02_input = input(' Chi ({0:.{1}f})? '.format(g_u02,PRE))
    if len(g_u02_input) != 0:
        g_u02 = float(g_u02_input)
    g_u03_input = input(' Phi ({0:.{1}f})? '.format(g_u03,PRE))
    if len(g_u03_input) != 0:
        g_u03 = float(g_u03_input)
    g_u04_input = input(' Mu ({0:.{1}f})? '.format(g_u04,PRE))
    if len(g_u04_input) != 0:
        g_u04 = float(g_u04_input)
    g_u05_input = input(' Gam ({0:.{1}f})? '.format(g_u05,PRE))
    if len(g_u05_input) != 0:
        g_u05 = float(g_u05_input)
    print ('\nEnter primary-reflection HKL coordinates:')
    g_h0_input = input(' H ({0:.{1}f})? '.format(g_h0,PRE))
    if len(g_h0_input) != 0:
        g_h0 = float(g_h0_input)
    g_k0_input = input(' K ({0:.{1}f})? '.format(g_k0,PRE))
    if len(g_k0_input) != 0:
        g_k0 = float(g_k0_input)
    g_l0_input = input(' L ({0:.{1}f})? '.format(g_l0,PRE))
    if len(g_l0_input) != 0:
        g_l0 = float(g_l0_input)
    UB()
    wh_refresh()

# Set H, K, L of secondary reflection
# Positions of secondary reflection are current positions
# Usage: or1() or or1(H,K,L)
def or1(*args):
    global g_u10, g_u11, g_u12, g_u13, g_u14, g_u15
    global g_h1, g_k1, g_l1
    global g_lambda1
    g_lambda1 = LAMBDA
    g_u10, g_u11, g_u12, g_u13, g_u14, g_u15 = (TTH, TH, CHI, PHI, MU, GAM)
    if len(args) == 0:
        print ('\nEnter secondary-reflection HKL coordinates:')
        g_h1_input = input(' H ({0:.{1}f})? '.format(g_h1,PRE))
        g_k1_input = input(' K ({0:.{1}f})? '.format(g_k1,PRE))
        g_l1_input = input(' L ({0:.{1}f})? '.format(g_l1,PRE))
        if len(g_h1_input) != 0:
            g_h1 = float(g_h1_input)
        if len(g_k1_input) != 0:
            g_k1 = float(g_k1_input)
        if len(g_l1_input) != 0:
            g_l1 = float(g_l1_input)
    elif len(args) == 3:
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        g_h1, g_k1, g_l1 = args
    else:
        print ('\nUsage:  or1()  or  or1(H,K,L)\n')
        return
    UB()
    wh_refresh()

# Set H, K, L and positions of secondary reflection
def setor1():
    global g_u10, g_u11, g_u12, g_u13, g_u14, g_u15
    global g_h1, g_k1, g_l1
    global g_lambda1
    g_lambda1 = LAMBDA
    print ('\nEnter secondary-reflection angles:')
    g_u10_input = input(' Two Theta ({0:.{1}f})? '.format(g_u10,PRE))
    if len(g_u10_input) != 0:
        g_u10 = float(g_u10_input)
    g_u11_input = input(' Theta ({0:.{1}f})? '.format(g_u11,PRE))
    if len(g_u11_input) != 0:
        g_u11 = float(g_u11_input)
    g_u12_input = input(' Chi ({0:.{1}f})? '.format(g_u12,PRE))
    if len(g_u12_input) != 0:
        g_u12 = float(g_u12_input)
    g_u13_input = input(' Phi ({0:.{1}f})? '.format(g_u13,PRE))
    if len(g_u13_input) != 0:
        g_u13 = float(g_u13_input)
    g_u14_input = input(' Mu ({0:.{1}f})? '.format(g_u14,PRE))
    if len(g_u14_input) != 0:
        g_u14 = float(g_u14_input)
    g_u15_input = input(' Gam ({0:.{1}f})? '.format(g_u15,PRE))
    if len(g_u15_input) != 0:
        g_u15 = float(g_u15_input)
    print ('\nEnter secondary-reflection HKL coordinates:')
    g_h1_input = input(' H ({0:.{1}f})? '.format(g_h1,PRE))
    if len(g_h1_input) != 0:
        g_h1 = float(g_h1_input)
    g_k1_input = input(' K ({0:.{1}f})? '.format(g_k1,PRE))
    if len(g_k1_input) != 0:
        g_k1 = float(g_k1_input)
    g_l1_input = input(' L ({0:.{1}f})? '.format(g_l1,PRE))
    if len(g_l1_input) != 0:
        g_l1 = float(g_l1_input)
    UB()
    wh_refresh()

# Check consistency at or0 and or1
wdesc ('or_check()','Check consistency of present or0 and or1 values')
def or_check():
    # Record original conditions
    global FLAG_WH
    o_FLAG_WH = FLAG_WH
    o_tth, o_th, o_chi, o_phi, o_mu, o_gam = (TTH, TH, CHI, PHI, MU, GAM)
    FLAG_WH = False
    # Deviation of H,K,L at or0
    mv(tth=g_u00,th=g_u01,chi=g_u02,phi=g_u03,mu=g_u04,gam=g_u05)
    dH0 = H - g_h0
    dK0 = K - g_k0
    dL0 = L - g_l0
    # Deviation of H,K,L at or1
    mv(tth=g_u10,th=g_u11,chi=g_u12,phi=g_u13,mu=g_u14,gam=g_u15)
    dH1 = H - g_h1
    dK1 = K - g_k1
    dL1 = L - g_l1
    print ('Check self consistency {} '.format(g_sample),end='')
    print (' with ({0:.{3}f}, {1:.{3}f}, {2:.{3}f} '.format(g_aa,g_bb,g_cc,PRE),end='')
    print (' {1:.{3}f}, {2:.{3}f}, {2:.{3}f})'.format(g_al,g_be,g_ga,2))
    print ('   At or0 ({0:.{6}f}, {1:.{6}f}, {2:.{6}f}): dH = {3:.{6}f}, dK = {4:.{6}f}, dL = {5:.{6}f}'.format(g_h0,g_k0,g_l0,dH0,dK0,dL0,PRE))
    print ('   At or1 ({0:.{6}f}, {1:.{6}f}, {2:.{6}f}): dH = {3:.{6}f}, dK = {4:.{6}f}, dL = {5:.{6}f}'.format(g_h1,g_k1,g_l1,dH1,dK1,dL1,PRE))
    print ('')
    # Return to original conditions
    mv(tth=o_tth,th=o_th,chi=o_chi,phi=o_phi,mu=o_mu,gam=o_gam)
    FLAG_WH = o_FLAG_WH

# Swap primary reflection and secondary reflection
wdesc('or_swap()','Swaps primary and secondary reference vectors/angles')
def or_swap():
    global g_u00, g_u01, g_u02, g_u03, g_u04, g_u05
    global g_h0, g_k0, g_l0
    global g_lambda0
    global g_u10, g_u11, g_u12, g_u13, g_u14, g_u15
    global g_h1, g_k1, g_l1
    global g_lambda1
    or0 = (g_lambda0, g_u00, g_u01, g_u02, g_u03, g_u04, g_u05, g_h0, g_k0, g_l0)
    or1 = (g_lambda1, g_u10, g_u11, g_u12, g_u13, g_u14, g_u15, g_h1, g_k1, g_l1)
    g_lambda0, g_u00, g_u01, g_u02, g_u03, g_u04, g_u05, g_h0, g_k0, g_l0 = or1
    g_lambda1, g_u10, g_u11, g_u12, g_u13, g_u14, g_u15, g_h1, g_k1, g_l1 = or0
    UB()
    wh_refresh()

# Redirecting to setfrozen
def setmode():
    print ('Redirecting to setfrozen')
    setfrozen()

# Print H, K, L, ALPHA, BETA, AZIMUTH, SA, Q
def wh_refresh():
    global OMEGA, SA, ABSQ, H, K, L, AZIMUTH, ALPHA, BETA
    OMEGA = TH - TTH/2
    # Scattering angle
    SA = 2*scbasic.thetaD_angle_2(TTH/2,MU,GAM)
    # Scattering vector
    ABSQ = scbasic.Q_length(LAMBDA,SA/2)
    # Azimuth reference vector
    az = np.array([[g_haz],[g_kaz],[g_laz]], dtype=float)
    hb = scbasic.CheckBr(LAMBDA,M_U,M_B,TTH/2,OMEGA,CHI,PHI,MU,GAM)
    H = hb[0][0]
    K = hb[1][0]
    L = hb[2][0]
    # When azimuth reference vector is surface normal, ALPHA and BETA are incident and exit angles
    AZIMUTH, ALPHA, BETA = scbasic.CheckPsiAlphainBetaout(M_U,M_B,az,TTH/2,OMEGA,CHI,PHI,MU,GAM)

# Show current positions
wdesc('wh()  (==wa())','Prints present HKL, etc.  (SA==Scattering Angle)')
def wh():
    print('')
    print ('H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(H,K,L,PRE))
    print ('|Q| = {0:.3f} nm-1  SA = {1:.{3}f} deg  at  LAMBDA = {2:.{4}f} A'.format(ABSQ*10,SA,LAMBDA,PRE,PRE+2))
    print ('AZ = ({0}, {1}, {2})  AZIMUTH = {3:.{6}f} deg  ALPHA = {4:.{6}f}  BETA = {5:.{6}f}'.format(g_haz,g_kaz,g_laz,AZIMUTH,ALPHA,BETA,PRE))
    print ('Omega = th-tth/2 = {0:.{1}f}'.format(OMEGA,PRE))
    print ('')
    strfmt = ('{:>'+str(PRE+6)+'}')*6
    strprt = ('tth','th','chi','phi','mu','gam')
    print (strfmt.format(*strprt))
    posfmt = ('{:>'+str(PRE+6)+'.'+str(PRE)+'f}')*6
    posprt = (TTH,TH,CHI,PHI,MU,GAM)
    print (posfmt.format(*posprt))
    print ('')

# Redirecting to wh()
def wa():
    print ('Redirecting to wh')
    wh()

# Start printing positions after mv() or br()
wdesc('wh_on(), wh_off()','turns on/off wh() call at end of br() and mv() calls'  )
def wh_on():
    global FLAG_WH
    FLAG_WH = True
    print ('')
    print ('Will print positions after mv() or br()')
    print ('')
# End printing positions after mv() or br()
def wh_off():
    global FLAG_WH
    FLAG_WH = False
    print ('')
    print ('End printing positions after mv() or br()')
    print ('')

# Change motors positions
# Usage: mv(mu=?,gam=?,tth=?,th=?,chi=?,phi=?)
wdesc('mv(tth=?,th=?,chi=?,phi=?,mu=?,gam=?)','Change angles')
def mv(**args):
    if len(args) == 0:
        print ('\nUsage: mv(tth=?,th=?,chi=?,phi=?,mu=?,gam=?)\n')
        return
    global TTH, TH, CHI, PHI, MU, GAM
    for key in args.keys():
        if (key in mnemonics) == False:
            print ('Invalid motor mnemonic for mv: {0}\n'.format(key))
            print ('Valid motor mnemonic for mv: tth, th, chi, phi, mu, gam')
            return
    for value in args.values():
        if (type(value) in [int,float]) == False:
            print ('Invalid motor position for mv: {0}\n'.format(value))
            return
    if 'tth' in args.keys():
        TTH = args['tth']
    if 'th' in args.keys():
        TH = args['th']
    if 'chi' in args.keys():
        CHI = args['chi']
    if 'phi' in args.keys():
        PHI = args['phi']
    if 'mu' in args.keys():
        MU = args['mu']
    if 'gam' in args.keys():
        GAM = args['gam']
    # Get (H, K, L, SA, OMEGA, ALPHA, BETA, AZIMUTH, ABSQ) at new positions
    wh_refresh()
    if FLAG_WH == True:
        wh()


# Set limits of positions
# Usage:  setlm()  or  setlm(ltth=?,utth=?,lth=?,uth=?,lchi=?,uchi=?,lphi=?,uphi=?,lmu=?,umu=?,lgam=?,ugam=?,lalpha=?,ualpha=?,lbeta=?,ubeta=?)
wdesc("setlm() , showlm()",'Sets/Shows limits.  Useful selecting amoung multiple solutions for (H,K,L) ')
def setlm(**args):
    global L_TTH, U_TTH, L_TH, U_TH, L_CHI, U_CHI, L_PHI, U_PHI, L_MU, U_MU, L_GAM, U_GAM, L_ALPHA, U_ALPHA, L_BETA, U_BETA
    if len(args) == 0:
        print ('\nSet limit of positions:')
        print ('')
        L_TTH_input = input(' Lower limit of tth ({0:.{1}f})? '.format(L_TTH,PRE))
        if len(L_TTH_input) != 0:
            L_TTH = float(L_TTH_input)
        U_TTH_input = input(' Upper limit of tth ({0:.{1}f})? '.format(U_TTH,PRE))
        if len(U_TTH_input) != 0:
            U_TTH = float(U_TTH_input)
        print ('')
        L_TH_input = input(' Lower limit of th ({0:.{1}f})? '.format(L_TH,PRE))
        if len(L_TH_input) != 0:
            L_TH = float(L_TH_input)
        U_TH_input = input(' Upper limit of th ({0:.{1}f})? '.format(U_TH,PRE))
        if len(U_TH_input) != 0:
            U_TH = float(U_TH_input)
        print ('')
        L_CHI_input = input(' Lower limit of chi ({0:.{1}f})? '.format(L_CHI,PRE))
        if len(L_CHI_input) != 0:
            L_CHI = float(L_CHI_input)
        U_CHI_input = input(' Upper limit of chi ({0:.{1}f})? '.format(U_CHI,PRE))
        if len(U_CHI_input) != 0:
            U_CHI = float(U_CHI_input)
        print ('')
        L_PHI_input = input(' Lower limit of phi ({0:.{1}f})? '.format(L_PHI,PRE))
        if len(L_PHI_input) != 0:
            L_PHI = float(L_PHI_input)
        U_PHI_input = input(' Upper limit of phi ({0:.{1}f})? '.format(U_PHI,PRE))
        if len(U_PHI_input) != 0:
            U_PHI = float(U_PHI_input)
        print ('')
        L_MU_input = input(' Lower limit of mu ({0:.{1}f})? '.format(L_MU,PRE))
        if len(L_MU_input) != 0:
            L_MU = float(L_MU_input)
        U_MU_input = input(' Upper limit of mu ({0:.{1}f})? '.format(U_MU,PRE))
        if len(U_MU_input) != 0:
            U_MU = float(U_MU_input)
        print ('')
        L_GAM_input = input(' Lower limit of gam ({0:.{1}f})? '.format(L_GAM,PRE))
        if len(L_GAM_input) != 0:
            L_GAM = float(L_GAM_input)
        U_GAM_input = input(' Upper limit of gam ({0:.{1}f})? '.format(U_GAM,PRE))
        if len(U_GAM_input) != 0:
            U_GAM = float(U_GAM_input)
        print ('')
        L_ALPHA_input = input(' Lower limit of alpha ({0:.{1}f})? '.format(L_ALPHA,PRE))
        if len(L_ALPHA_input) != 0:
            L_ALPHA = float(L_ALPHA_input)
        U_ALPHA_input = input(' Upper limit of alpha ({0:.{1}f})? '.format(U_ALPHA,PRE))
        if len(U_ALPHA_input) != 0:
            U_ALPHA = float(U_ALPHA_input)
        print ('')
        L_BETA_input = input(' Lower limit of beta ({0:.{1}f})? '.format(L_BETA,PRE))
        if len(L_BETA_input) != 0:
            L_BETA = float(L_BETA_input)
        U_BETA_input = input(' Upper limit of beta ({0:.{1}f})? '.format(U_BETA,PRE))
        if len(U_BETA_input) != 0:
            U_BETA = float(U_BETA_input)
        print ('')
    else:
        lmkeys = ['ltth','utth','lth','uth','lchi','uchi','lphi','uphi','lmu','umu','lgam','ugam','lalpha','ualpha','lbeta','ubeta']
        for key in args.keys():
            if (key in lmkeys) == False:
                print ('Invalid keyword for setlm: {0}'.format(key))
                print ('Valid keyword for setlm: l(u)tth, l(u)th, l(u)chi, l(u)phi, l(u)mu, l(u)gam, l(u)alpha, l(u)beta')
                return
        for value in args.values():
            if (type(value) in [int,float]) == False:
                print ('Invalid value for setlm: {0}\n'.format(value))
                return
        if 'ltth' in args.keys():
            L_TTH = args['ltth']
        if 'utth' in args.keys():
            U_TTH = args['utth']
        if 'lth' in args.keys():
            L_TH = args['lth']
        if 'uth' in args.keys():
            U_TH = args['uth']
        if 'lchi' in args.keys():
            L_CHI = args['lchi']
        if 'uchi' in args.keys():
            U_CHI = args['uchi']
        if 'lphi' in args.keys():
            L_PHI = args['lphi']
        if 'uphi' in args.keys():
            U_PHI = args['uphi']
        if 'lmu' in args.keys():
            L_MU = args['lmu']
        if 'umu' in args.keys():
            U_MU = args['umu']
        if 'lgam' in args.keys():
            L_GAM = args['lgam']
        if 'ugam' in args.keys():
            U_GAM = args['ugam']
        if 'lalpha' in args.keys():
            L_ALPHA = args['lalpha']
        if 'lalpha' in args.keys():
            U_ALPHA = args['ualpha']
        if 'lbeta' in args.keys():
            L_BETA = args['lbeta']
        if 'ubeta' in args.keys():
            U_BETA = args['ubeta']
        showmainlm()

def showlm():
    print("\n  Motor limits are now",end='')
    print("    tth: {} {} ".format(L_TTH,U_TTH), end='')
    print("     th: {} {} ".format(L_TH,U_TH), end='')
    print("    chi: {} {} ".format(L_CHI,U_CHI), end='')
    print("    phi: {} {} ".format(L_PHI,U_PHI), end='')
    print("     mu: {} {} ".format(L_MU,U_MU), end='')
    print("    gam: {} {} ".format(L_GAM,U_GAM), end='')
    print("  alpha: {} {} ".format(L_ALPHA,U_ALPHA), end='')
    print("   beta: {} {} ".format(L_BETA,U_BETA), end='')
    print("\n")

def showmainlm():
    print("\n  Main motor limits are now",end='')
    print("    tth: {} {} ".format(L_TTH,U_TTH), end='')
    print("     th: {} {} ".format(L_TH,U_TH), end='')
    print("    chi: {} {} ".format(L_CHI,U_CHI), end='')
    print("    phi: {} {} ".format(L_PHI,U_PHI), end='')
    print("\n")

# Set all limits as -180 to 180
wdesc("setlm_clear()",'Sets default (+-180 deg) limits on all angles')
def setlm_clear():
    global L_TTH, U_TTH, L_TH, U_TH, L_CHI, U_CHI, L_PHI, U_PHI, L_MU, U_MU, L_GAM, U_GAM, L_ALPHA, U_ALPHA, L_BETA, U_BETA
    L_TTH = -180.0; U_TTH = 180.0
    L_TH = -180.0; U_TH = 180.0
    L_CHI = -180.0; U_CHI = 180.0 
    L_PHI = -180.0; U_PHI = 180.0
    L_MU = -180.0; U_MU = 180.0
    L_GAM = -180.0; U_GAM = 180.0
    L_ALPHA = -180.0; U_ALPHA = 180.0 
    L_BETA = -180.0; U_BETA = 180.0
    print ('')
    print ('Now all limits are set as -180 to 180')
    print ('')

# Calculate positions in specified frozen
# Show the 1st set of positions in preset limits
# Usage:  ca(H,K,L)
wdesc("ca(H,K,L)","Calculate angles for a given (H,K,L)")
def ca(*args):
    if len(args) != 3:
        print ('Usage:  ca(H,K,L)')
        return
    for value in args:
        if (type(value) in [int,float]) == False:
            print ('\nInvalid argument: {0}\n'.format(value))
            return
    caH, caK, caL = args
    flag, pos = ca_s(caH, caK, caL)
    if flag == False:
        dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
        dic_pos = {0:F_TTH, 1:F_TH, 2:F_CHI, 3:F_PHI, 4:F_MU, 5:F_GAM, 6:F_OMEGA, 7:F_AZIMUTH, 8:F_ALPHA, 9:F_BETA}
        g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(g_frozen)]
        print ('Error: Impossible reflection within current limits for frozen: ', end='')
        freprt = (g_frozen,dic_ang[g_frozen_1],dic_pos[g_frozen_1],dic_ang[g_frozen_2],dic_pos[g_frozen_2],dic_ang[g_frozen_3],dic_pos[g_frozen_3])
        print ('{1}={2:.{7}f} {3}={4:.{7}f} {5}={6:.{7}f}'.format(*freprt,PRE))
        return
    else:
        caTTH, caTH, caCHI, caPHI, caMU, caGAM, caSA, caOMEGA, caAZIMUTH, caALPHA, caBETA = pos[0]
        caABSQ = scbasic.Q_length(LAMBDA,abs(caSA/2))
        print('')
        print ('H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(caH,caK,caL,PRE))
        print ('|Q| = {0:.3f} nm-1  SA = {1:.{3}f} deg  at  LAMBDA = {2:.{4}f} A'.format(caABSQ*10,caSA,LAMBDA,PRE,PRE+2))
        print ('AZ = ({0}, {1}, {2})  AZIMUTH = {3:.{6}f} deg  ALPHA = {4:.{6}f}  BETA = {5:.{6}f}'.format(g_haz,g_kaz,g_laz,caAZIMUTH,caALPHA,caBETA,PRE))
        print ('Omega = th-tth/2 = {0:.{1}f}'.format(caOMEGA,PRE))
        print ('')
        strfmt = ('{:>'+str(PRE+6)+'}')*6
        strprt = ('tth','th','chi','phi','mu','gam')
        print (strfmt.format(*strprt))
        posfmt = ('{:>'+str(PRE+6)+'.'+str(PRE)+'f}')*6
        posprt = (caTTH,caTH,caCHI,caPHI,caMU,caGAM)
        print (posfmt.format(*posprt))
        print ('')
        print ('Command (sixcircle):  ', end='')
        print ('mv (tth={0:.{6}f}, th={1:.{6}f}, chi={2:.{6}f}, phi={3:.{6}f}, mu={4:.{6}f}, gam={5:.{6}f})'.format(*posprt,PRE))
        print ('Command (BL43LXU):    ', end='')
        print ('mv tth {0:.{6}f} th {1:.{6}f} chi {2:.{6}f} phi {3:.{6}f}'.format(*posprt,PRE))
        print ('')

# Move to H, K, L    Using first set of positions from ca_s()
# Usage: br(H,K,L)
wdesc("br(H,K,L)","Move to given Q=(H,K,L)")
def br(*args):
    if len(args) != 3:
        print ('Usage:  br(H,K,L)')
        return
    for value in args:
        if (type(value) in [int,float]) == False:
            print ('\nInvalid argument: {0}\n'.format(value))
            return
    goH, goK, goL = args
    flag, pos = ca_s(goH, goK, goL)
    dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
    dic_pos = {0:F_TTH, 1:F_TH, 2:F_CHI, 3:F_PHI, 4:F_MU, 5:F_GAM, 6:F_OMEGA, 7:F_AZIMUTH, 8:F_ALPHA, 9:F_BETA}
    g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(g_frozen)]
    if flag == False:
        print ('Error: Impossible reflection within current limits for frozen {0}: ', end='')
        freprt = (g_frozen,dic_ang[g_frozen_1],dic_pos[g_frozen_1],dic_ang[g_frozen_2],dic_pos[g_frozen_2],dic_ang[g_frozen_3],dic_pos[g_frozen_3])
        print ('freeze {1}={2:.{7}f} {3}={4:.{7}f} {5}={6:.{7}f}'.format(*freprt,PRE))
        return
    goprt = (goH,goK,goL,dic_ang[g_frozen_1],dic_ang[g_frozen_2],dic_ang[g_frozen_3],dic_pos[g_frozen_1],dic_pos[g_frozen_2],dic_pos[g_frozen_3])
    print ('Moving to ({0}, {1}, {2}) with ({3}, {4}, {5}) frozen at ({6:.{9}f}, {7:.{9}f}, {8:.{9}f})'.format(*goprt,PRE))
    go_tth,go_th,go_chi,go_phi,go_mu,go_gam,go_sa,go_omega,go_azimuth,go_alpha,go_beta = pos[0]
    mv(tth=go_tth,th=go_th,chi=go_chi,phi=go_phi,mu=go_mu,gam=go_gam)

# Calculate positions in specified frozen
# _s: silent, calculation in background
# return to: flag, pos
# flag: True (False) when calculation is successful (failed: no solution)
# pos: N set of positions, [[tth1, th1, chi1, phi1, mu1, gam1, sa1, omega1, azimuth1, alpha1, beta1],[tth2,...],...,[tthN,...]]
def ca_s(caH, caK, caL):
    # Azimuth reference vector
    az = np.array([[g_haz],[g_kaz],[g_laz]], dtype=float)
    # Target H, K, L
    ca_hb = np.array([[caH],[caK],[caL]], dtype=float)
    # Pre-specified frozen and positions of frozen angles
    dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
    dic_pos = {0:F_TTH, 1:F_TH, 2:F_CHI, 3:F_PHI, 4:F_MU, 5:F_GAM, 6:F_OMEGA, 7:F_AZIMUTH, 8:F_ALPHA, 9:F_BETA}
    g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(g_frozen)]
    # Prepare for calculation. Put in positions of specified frozen angles. For unfrozen angles to be calculated: 'x'
    dic = {'tth':'x','th':'x','chi':'x','phi':'x','mu':'x','gam':'x','omega':'x','azimuth':'x','alpha':'x','beta':'x','theta':'x'}
    dic[dic_ang[g_frozen_1]] = dic_pos[g_frozen_1]
    dic[dic_ang[g_frozen_2]] = dic_pos[g_frozen_2]
    dic[dic_ang[g_frozen_3]] = dic_pos[g_frozen_3]
    # If freeze tth, use theta = tth/2 in calculation. Note theta != th
    g_frozen_str = str(g_frozen)
    if g_frozen_str.count('0') == 1:
        dic['theta'] = dic['tth']/2
    flag = False
    pos = []
    # In general: neither azimuth, alpha, beta is frozen
    if g_frozen_str.count('7') + g_frozen_str.count('8') + g_frozen_str.count('9') == 0:
        flagp, angles = scbasic.SC_angles(LAMBDA,M_U,M_B,ca_hb,dic['th'],dic['theta'],dic['omega'],dic['chi'],dic['phi'],dic['mu'],dic['gam'])
    else:
        # Azimuth frozen
        if g_frozen_str.count('7') == 1:
            flagp, angles = scbasic.SC_angles_fix_psi(LAMBDA,M_U,M_B,ca_hb,az,dic['azimuth'],dic['th'],dic['theta'],dic['omega'],dic['chi'],dic['phi'],dic['mu'],dic['gam'])
        # Alpha frozen or Beta frozen
        else:
            flagp, angles = scbasic.SC_angles_fix_ab(LAMBDA,M_U,M_B,ca_hb,az,dic['alpha'],dic['beta'],dic['th'],dic['theta'],dic['omega'],dic['chi'],dic['phi'],dic['mu'],dic['gam'])
    if flagp == True:
        # Current angles: [[theta,omega,chi,phi,mu,gam],...]
        # Sort out angles to [[tth,th,chi,phi,mu,gam],...]
        N, angles = scbasic.motors(angles)
        for element in angles:
            caTTH, caTH, caCHI, caPHI, caMU, caGAM = element
            caSA = 2*scbasic.thetaD_angle_2(caTTH/2,caMU,caGAM)
            caOMEGA = caTH - caTTH/2
            # Calculate corresponding azimuth, alpha, beta
            caAZIMUTH, caALPHA, caBETA = scbasic.CheckPsiAlphainBetaout(M_U,M_B,az,caTTH/2,caOMEGA,caCHI,caPHI,caMU,caGAM)
            if L_TTH <= caTTH < U_TTH and L_TH <= caTH < U_TH and L_CHI <= caCHI < U_CHI and L_PHI <= caPHI < U_PHI:
                if L_MU <= caMU < U_MU and L_GAM <= caGAM < U_GAM:
                    if L_ALPHA <= caALPHA < U_ALPHA and L_BETA <= caBETA < U_BETA:
                        pos.append([caTTH,caTH,caCHI,caPHI,caMU,caGAM,caSA,caOMEGA,caAZIMUTH,caALPHA,caBETA])
    if len(pos) > 0:
        flag = True
    return flag, pos

# Calculate positions in specified frozen
# _a: all, show all sets of positions in preset limits
# Usage:  ca_a(H,K,L)
wdesc("ca_a(H,K,L), ca_s(H,K,L) ","Calculates (_s=silently) angles for (H,K,L) within limits from setlm ")
def ca_a(*args):
    if len(args) != 3:
        print ('Usage:  ca_a(H,K,L)')
        return
    for value in args:
        if (type(value) in [int,float]) == False:
            print ('\nInvalid argument: {0}\n'.format(value))
            return
    caH, caK, caL = args
    flag, pos = ca_s(caH, caK, caL)
    if flag == False:
        dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
        dic_pos = {0:F_TTH, 1:F_TH, 2:F_CHI, 3:F_PHI, 4:F_MU, 5:F_GAM, 6:F_OMEGA, 7:F_AZIMUTH, 8:F_ALPHA, 9:F_BETA}
        g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(g_frozen)]
        print ('Error: Impossible reflection within current limits for frozen {0}: ', end='')
        freprt = (g_frozen,dic_ang[g_frozen_1],dic_pos[g_frozen_1],dic_ang[g_frozen_2],dic_pos[g_frozen_2],dic_ang[g_frozen_3],dic_pos[g_frozen_3])
        print ('freeze {1}={2:.{7}f} {3}={4:.{7}f} {5}={6:.{7}f}'.format(*freprt,PRE))
        return
    else:
        caSA = pos[0][6]
        caABSQ = scbasic.Q_length(LAMBDA,abs(caSA/2))
        print ('')
        print ('Calculated Positions:')
        print ('')
        print ('H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(caH,caK,caL,PRE))
        print ('|Q| = {0:.3f} nm-1  SA = {1:.{3}f} deg  at  LAMBDA = {2:.{4}f} A'.format(caABSQ*10,caSA,LAMBDA,PRE,PRE+2))
        print ('')
        strfmt = '    '+('{:>'+str(PRE+6)+'}')*11
        strprt = ('tth','th','chi','phi','mu','gam','sa','omega','azimuth','alpha','beta')
        print (strfmt.format(*strprt))
        posfmt = '{:>4}'+('{:>'+str(PRE+6)+'.'+str(PRE)+'f}')*11
        for i in range(0,len(pos)):
            print (posfmt.format(i,*pos[i],PRE))
        print ('')
        print ('Command (sixcircle):')
        for i in range(0,len(pos)):
            posprt = [pos[i][0],pos[i][1],pos[i][2],pos[i][3],pos[i][4],pos[i][5]]
            print ('{0:>4}  mv (tth={1:.{7}f}, th={2:.{7}f}, chi={3:.{7}f}, phi={4:.{7}f}, mu={5:.{7}f}, gam={6:.{7}f})'.format(i,*posprt,PRE))        
        print ('')


# Check where the limits of alpha and beta are for H, K, L
# _s: silent, calculation in background
# Usage:  flag, min, max = wmab_s(H, K, L)
def wmab_s(caH, caK, caL):
    # Azimuth reference vector
    az = np.array([[g_haz],[g_kaz],[g_laz]], dtype=float)
    # Target H, K, L
    ca_hb = np.array([[caH],[caK],[caL]], dtype=float)
    flag, minab, maxab = scbasic.CheckRangeAlphainBetaout(LAMBDA,M_B,ca_hb,az)
    return flag, minab, maxab

# Check where the limits of alpha and beta are for H, K, L
# Usage:  wmab(H, K, L)
wdesc("wmab(H,K,L), wmab_s(H,K,L)","Check the limits on alpha and beta for given (H,K,L)")
def wmab(*args):
    if len(args) != 3:
        print ('Usage:  wmab(H,K,L)')
        return
    for value in args:
        if (type(value) in [int,float]) == False:
            print ('\nInvalid argument: {0}\n'.format(value))
            return
    caH, caK, caL = args
    flag, minab, maxab = wmab_s(caH, caK, caL)
    if flag == False:
        print ('')
        print ('Error: Impossible reflection.')
        return
    print ('')
    print ('Limits of ALPHA and BETA for {0}, {1}, {2}'.format(caH, caK, caL))
    print ('')
    print ('{0:>10}{1:>10}'.format('Min','Max'))
    print ('{0:>10.{2}f}{1:>10.{2}f}'.format(minab,maxab,PRE))
    print ('')

# Set output precision
# Usage:  setprecision()  or setprecision(n)
# Default n = 4
wdesc("setprecision(n)","Sets decimal precision for output (default n=4)")
def setprecision(*args):
    global PRE
    if len(args) == 0:
        PRE_input =  input('\nOutput precision ({0})? '.format(PRE))
        if len(PRE_input) != 0:
            PRE = int(PRE_input)
    elif len(args) == 1:
        if type(args[0]) != int:
            print ('\nInvalid argument: {0}\n'.format(args[0]))
            return
        if args[0] < 0:
            print ('\nInvalid argument: {0}\n'.format(args[0]))
            return          
        PRE = int(args[0])
    else:
        print ('\nUsage:  setprecision()  or  setprecision(n)\n')
        return
    print ('')
    print ('Output precision set to {0}'.format(PRE))
    print ('')

ini()
