# sixcircle: code for six-circle diffractometer
# Required additional documents: 1- scbasic.py, 2- ini.conf

# Developed with Python 3.7.3, numpy 1.16.4, scipy 1.3.0
# Materials Dynamics Laboratory, RIKEN SPring-8 Center
# Ver 1.51, December 2020

# Authors: Wenyang Zhao & Alfred Q. R. Baron
# Contact: baron@spring8.or.jp

# If this code is used independently of work at BL43, please reference the writeup ("Open-source Python software for six-circle diffraction with an inelastic x-ray scattering (IXS) spectrometer" by Wenyang ZHAO and Alfred Q.R. BARON, unpublished, available at https://beamline.harima.riken.jp/bl43lxu)

# In typical use related to experimental work at BL43LXU, it is enough to reference a BL publication.  However, if this code is used extensively, then please include a separate reference as above.

# 04/23/2025: This code was refactored to make it run in a bluesky environment -- Juan Marulanda

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
import scbasic


class SixCircle():
    
    def __init__(self):
        self.runquiet = False
        
        self.ini()

    def __doc__(self):
        print('Function Definitions from sixcicle.py:')
        self.wdesc(self.ini)
        self.wdesc(self.load)
        self.wdesc(self.save)
        self.wdesc(self.UB)
        self.wdesc(self.pa)
        self.wdesc(self.setlambda)
        self.wdesc(self.setlat)
        self.wdesc(self.setfrozen)
        self.wdesc(self.freeze)
        self.wdesc(self.setaz)
        self.wdesc(self.or0)
        self.wdesc(self.setor0)
        self.wdesc(self.or_check)
        self.wdesc(self.or_swap)
        self.wdesc(self.wh)
        self.wdesc(self.wh_on)
        self.wdesc(self.mv)
        self.wdesc(self.setlm)
        self.wdesc(self.setlm_clear)
        self.wdesc(self.ca)
        self.wdesc(self.br)
        self.wdesc(self.ca_a)
        self.wdesc(self.wmab)
        self.wdesc(self.setprecision)
    
    def wdesc(self, fname):
        print(f'{fname.__name__}',end='')
        print(f'{fname.__doc__}')
    
    def runquiet_on(self):
        self.runquiet = True
    
    def runquiet_off(self):
        self.runquiet = False

    # Initialization
    def ini(self):
        """() \t\t\t\t\t\t\t\t\t\t SixCircle initialization"""
        
        print('\nInitilizing sixcircle setup.')
        # Valid mnemonics
        self.mnemonics = ['tth','th','chi','phi','mu','gam']
        # Positions: all zero at startup
        self.TTH=0.0; self.TH=0.0; self.CHI=0.0; self.PHI=0.0; self.MU=0.0; self.GAM=0.0
        # SA: scattering angle; ABSQ: length of scattering vector
        self.SA = 0.0; self.ABSQ = 0.0
        # Omega: difference from TTH/2 to TH
        self.OMEGA = self.TH - self.TTH/2
        # Frozen of six-circle calculation. Default '456' at startup (freeze mu, gam, omega)
        self.g_frozen = '456'
        # Frozen positions: all zero at startup
        self.F_TTH=0.0; self.F_TH=0.0; self.F_CHI=0.0; self.F_PHI=0.0; self.F_MU=0.0
        self.F_GAM=0.0; self.F_OMEGA=0.0; self.F_AZIMUTH=0.0; self.F_ALPHA=0.0; self.F_BETA=0.0
        # Output precision, default PRE=3, (PRE+2) for wavelength
        self.PRE = 3
        # Flag of using wh() after mv() or br(), default False
        self.FLAG_WH = True
        # Load default configuration file at startup
        if path.isfile('./conf/sixcircle_last_UB') :
            self.load('./conf/sixcircle_last_UB')
        else :
            self.load('./conf/ini.conf')
    
    # Load a configuration file
    def load(self, filepath):
        """('filepath') \t\t\t\t\t\t\t Load a configuration file. Loads Wavelength and orientation matrix"""
        try:
            with open(filepath, 'r') as f:
                content = f.readlines()
            if not self.runquiet : print ('Reading configuration from {0}'.format(filepath))
        except:
            print ('\nError reading configuration file {0}'.format(filepath))
            return
        
        dic = {}
        for line in content:
            if line.startswith('GLOBAL') == True:
                gvar = line.split()[1]
                gvarvalue = line.split()[2]
                dic.update({gvar:gvarvalue})

        # Sample description
        self.g_sample = dic['g_sample']
        # Azimuth reference vector H, K, L
        self.g_haz = float(dic['g_haz'])
        self.g_kaz = float(dic['g_kaz'])
        self.g_laz = float(dic['g_laz'])
        # Lattice parameters: a, b, c, alpha, beta, gam
        self.g_aa = float(dic['g_aa'])
        self.g_bb = float(dic['g_bb'])
        self.g_cc = float(dic['g_cc'])
        self.g_al = float(dic['g_al'])
        self.g_be = float(dic['g_be'])
        self.g_ga = float(dic['g_ga'])
        # Primary reflection: H, K, L
        self.g_h0 = float(dic['g_h0'])
        self.g_k0 = float(dic['g_k0'])
        self.g_l0 = float(dic['g_l0'])
        # Primary reflection: positions of tth, th, chi, phi, mu, gam
        self.g_u00 = float(dic['g_u00'])
        self.g_u01 = float(dic['g_u01'])
        self.g_u02 = float(dic['g_u02'])
        self.g_u03 = float(dic['g_u03'])
        self.g_u04 = float(dic['g_u04'])
        self.g_u05 = float(dic['g_u05'])
        # Secondary reflections: H, K, L
        self.g_h1 = float(dic['g_h1'])
        self.g_k1 = float(dic['g_k1'])
        self.g_l1 = float(dic['g_l1'])
        # Secondary reflection: positions of tth, th, chi, phi, mu, gam
        self.g_u10 = float(dic['g_u10'])
        self.g_u11 = float(dic['g_u11'])
        self.g_u12 = float(dic['g_u12'])
        self.g_u13 = float(dic['g_u13'])
        self.g_u14 = float(dic['g_u14'])
        self.g_u15 = float(dic['g_u15'])
        # Wavelength in finding reference reflections
        self.g_lambda0 = float(dic['g_lambda0'])
        self.g_lambda1 = float(dic['g_lambda1'])
        # Limit of positions
        self.L_TTH = float(dic['L_TTH'])
        self.U_TTH = float(dic['U_TTH'])
        self.L_TH = float(dic['L_TH'])
        self.U_TH = float(dic['U_TH'])
        self.L_CHI = float(dic['L_CHI'])
        self.U_CHI = float(dic['U_CHI'])
        self.L_PHI = float(dic['L_PHI'])
        self.U_PHI = float(dic['U_PHI'])
        self.L_MU = float(dic['L_MU'])
        self.U_MU = float(dic['U_MU'])
        self.L_GAM = float(dic['L_GAM'])
        self.U_GAM = float(dic['U_GAM'])
        self.L_ALPHA = float(dic['L_ALPHA'])
        self.U_ALPHA = float(dic['U_ALPHA'])
        self.L_BETA = float(dic['L_BETA'])
        self.U_BETA = float(dic['U_BETA'])
        # # Set wavelength in current calculation
        self.LAMBDA = self.g_lambda0
        # # Calculate matrix U and B
        self.UB()
        
    # Save a configuration file
    def save(self, filepath):
        """('filepath') \t\t\t\t\t\t\t Save a configuration file. Saves crystal parameters, Wavelength, and orientation matrix"""
        
        try:
            with open(filepath, 'w') as f:
                f.write('# Configuration file of sixcircle.\n')
                f.write('\n')
                f.write('# Sample description\n')
                f.write('GLOBAL g_sample {0}\n'.format(self.g_sample))
                f.write('\n')
                f.write('# Azimuthal reference H K L\n')
                dic = dict(g_haz=self.g_haz, g_kaz=self.g_kaz, g_laz=self.g_laz)
                for key in dic.keys():
                    f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],self.PRE))
                f.write('\n')
                f.write('# Lattice parameters\n')
                dic = dict(g_aa=self.g_aa, g_bb=self.g_bb, g_cc=self.g_cc, g_al=self.g_al, g_be=self.g_be, g_ga=self.g_ga)
                for key in dic.keys():
                    f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],self.PRE))
                f.write('\n')
                f.write('# Primary-reflection wavelength\n')
                f.write('GLOBAL g_lambda0 {0:.{1}f}\n'.format(self.g_lambda0,self.PRE+2))
                f.write('\n')
                f.write('# Primary-reflection HKL coordinates\n')
                dic = dict(g_h0=self.g_h0, g_k0=self.g_k0, g_l0=self.g_l0)
                for key in dic.keys():
                    f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],self.PRE))
                f.write('\n')
                f.write('# Primary-reflection angles: tth, th, chi, phi, mu, gam\n')
                dic = dict(g_u00=self.g_u00, g_u01=self.g_u01, g_u02=self.g_u02, g_u03=self.g_u03, g_u04=self.g_u04, g_u05=self.g_u05)
                for key in dic.keys():
                    f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],self.PRE))
                f.write('\n')
                f.write('# Secondary-reflection wavelength\n')
                f.write('GLOBAL g_lambda1 {0:.{1}f}\n'.format(self.g_lambda1,self.PRE+2))
                f.write('\n')
                f.write('# Secondary-reflection HKL coordinates\n')
                dic = dict(g_h1=self.g_h1, g_k1=self.g_k1, g_l1=self.g_l1)
                for key in dic.keys():
                    f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],self.PRE))
                f.write('\n')
                f.write('# Secondary-reflection angles: tth, th, chi, phi, mu, gam\n')
                dic = dict(g_u10=self.g_u10, g_u11=self.g_u11, g_u12=self.g_u12, g_u13=self.g_u13, g_u14=self.g_u14, g_u15=self.g_u15)
                for key in dic.keys():
                    f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic[key],self.PRE))
                f.write('\n')
                f.write('# Limit of positions\n')
                dic1 = dict(L_TTH=self.L_TTH, U_TTH=self.U_TTH, L_TH=self.L_TH, U_TH=self.U_TH, L_CHI=self.L_CHI, U_CHI=self.U_CHI, L_PHI=self.L_PHI, U_PHI=self.U_PHI)
                dic2 = dict(L_MU=self.L_MU, U_MU=self.U_MU, L_GAM=self.L_GAM, U_GAM=self.U_GAM, L_ALPHA=self.L_ALPHA, U_ALPHA=self.U_ALPHA, L_BETA=self.L_BETA, U_BETA=self.U_BETA)
                for key in dic1.keys():
                    f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic1[key],self.PRE))
                for key in dic2.keys():
                    f.write('GLOBAL {0} {1:.{2}f}\n'.format(key,dic2[key],self.PRE))
            if not self.runquiet : print ('\nWrote configuration to {0}'.format(filepath))
        except:
            print ('\nError in writing configuration file {0}'.format(filepath))

    def UB(self):
        """() \t\t\t\t\t\t\t\t\t\t Updates UB matrix"""
        
        # wavelength in finding reference reflections
        tupleB = scbasic.B_matrix(self.g_aa, self.g_bb, self.g_cc, self.g_al, self.g_be, self.g_ga)
        if (tupleB[0]) == False:
            print ('\nError! Invalid lattice parameters!')
            return
        # Lattice parameters in reciprocal space: a, b, c, alpha, beta, gam
        flag, self.M_B, self.g_aa_s, self.g_bb_s, self.g_cc_s, self.g_al_s, self.g_be_s, self.g_ga_s = tupleB
        # Primary and secondary reflections: unit vectors calcualted from positions
        u0phi = scbasic.uphi_vector(self.g_u00/2, self.g_u01-self.g_u00/2, self.g_u02, self.g_u03, self.g_u04, self.g_u05)
        u1phi = scbasic.uphi_vector(self.g_u10/2, self.g_u11-self.g_u10/2, self.g_u12, self.g_u13, self.g_u14, self.g_u15)
        # Primary and secondary reflections: unit vectors cacluated from H, K, L
        h0b = np.array([[self.g_h0],[self.g_k0],[self.g_l0]],dtype=float)
        h1b = np.array([[self.g_h1],[self.g_k1],[self.g_l1]],dtype=float)
        u0c = scbasic.uc_vector(self.M_B, h0b)
        u1c = scbasic.uc_vector(self.M_B, h1b)
        # Calculate orientation matrix U
        errcode, self.M_U = scbasic.U_matrix(u0phi, u1phi, u0c, u1c)
        if errcode == 0:
           # print ('\nUB recalculated for {} '.format(self.g_sample),end='')
           # print (' with ({0:.{3}f}, {1:.{3}f}, {2:.{3}f} '.format(self.g_aa,self.g_bb,self.g_cc,self.PRE),end='')
           # print (' {1:.{3}f}, {2:.{3}f}, {2:.{3}f})'.format(self.g_al,self.g_be,self.g_ga,2))
           # print ('     or0 = ({0:.{3}f}, {1:.{3}f}, {2:.{3}f})    and '.format(self.g_h0,self.g_k0,self.g_l0,self.PRE),end='')
           # print ('     or1 = ({0:.{3}f}, {1:.{3}f}, {2:.{3}f}) '.format(self.g_h1,self.g_k1,self.g_l1,self.PRE))
           ####################
           self.wh_refresh()
           self.or_check()
           self.runquiet_on()
           self.save('./conf/sixcircle_last_UB')
           self.runquiet_off()
        elif errcode == 1:
            print ('\nCannot find orientation matrix:  Reflections are parallel.\n')
        elif errcode == 2:
            print ('\nCannot find orientation matrix:  Reflections (by angles) are parallel.\n')
    
    def pa(self):
        """() \t\t\t\t\t\t\t\t\t\t Show parameters of orientation calculation. Display crystal structure and reference vectors"""
        
        print ('')
        print ('Primary Reflection (or0, at lambda {0:.{1}f}):'.format(self.g_lambda0,self.PRE+2))
        print ('           tth={0:.{6}f}, th={1:.{6}f}, chi={2:.{6}f}, phi={3:.{6}f}, mu={4:.{6}f}, gam={5:.{6}f}'.format(self.g_u00,self.g_u01,self.g_u02,self.g_u03,self.g_u04,self.g_u05,self.PRE))
        print ('           H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(self.g_h0,self.g_k0,self.g_l0,self.PRE))
        print ('')
        print ('Secondary Reflection (or1, at lambda {0:.{1}f}):'.format(self.g_lambda1,self.PRE+2))
        print ('           tth={0:.{6}f}, th={1:.{6}f}, chi={2:.{6}f}, phi={3:.{6}f}, mu={4:.{6}f}, gam={5:.{6}f}'.format(self.g_u10,self.g_u11,self.g_u12,self.g_u13,self.g_u14,self.g_u15,self.PRE))
        print ('           H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(self.g_h1,self.g_k1,self.g_l1,self.PRE))
        print ('')
        print ('Lattice Constants (lengths / angles):')
        print ('           real space = {0:.{6}f} {1:.{6}f} {2:.{6}f} / {3:.{6}f} {4:.{6}f} {5:.{6}f}'.format(self.g_aa,self.g_bb,self.g_cc,self.g_al,self.g_be,self.g_ga,self.PRE))
        print ('           reciprocal space = {0:.{6}f} {1:.{6}f} {2:.{6}f} / {3:.{6}f} {4:.{6}f} {5:.{6}f}'.format(self.g_aa_s,self.g_bb_s,self.g_cc_s,self.g_al_s,self.g_be_s,self.g_ga_s,self.PRE))
        print ('')
        print ('Azimuthal Reference:')
        print ('           H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(self.g_haz, self.g_kaz, self.g_laz,self.PRE))
        print ('')
        print ('           LAMBDA = {0:.{1}f}'.format(self.LAMBDA,self.PRE+2))
        print ('')

    # Set wavelength
    def setlambda(self, *args):
        """(wavelength) \t\t\t\t\t\t Sets LAMBDA in angstroms"""
        
        if len(args) == 0:
            LAMBDA_input =  input('\nWavelength / A ({0:.{1}f})? '.format(self.LAMBDA,self.PRE+2))
            if len(LAMBDA_input) != 0:
                self.LAMBDA = float(LAMBDA_input)
        elif len(args) == 1:
            if (type(args[0]) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(args[0]))
                return
            self.LAMBDA = float(args[0])
        else:
            print ('\nUsage:  setlambda()  or  setlambda(LAMBDA)\n')
            return
        print("  -> LAMBDA now %.8f A" %(self.LAMBDA))
        self.wh_refresh()
    
    # Usage: setlat() or setlat(a,b,c,alpha,beta,gam)
    def setlat(self, *args):
        """(), setlat(a,b,c,alpha,beta,gamma) \t Set lattice parameters. Set crystal parameters (A & deg) """
        
        if len(args) == 0:
            print ('\nEnter real space lattice parameters:')
            g_aa_input = input(' Lattice a ({0:.{1}f})? '.format(self.g_aa,self.PRE))
            if len(g_aa_input) != 0:
                self.g_aa = float(g_aa_input)
            g_bb_input = input(' Lattice b ({0:.{1}f})? '.format(self.g_bb,self.PRE))
            if len(g_bb_input) != 0:
                self.g_bb = float(g_bb_input)
            g_cc_input = input(' Lattice c ({0:.{1}f})? '.format(self.g_cc,self.PRE))
            if len(g_cc_input) != 0:
                self.g_cc = float(g_cc_input)
            g_al_input = input(' Lattice alpha ({0:.{1}f})? '.format(self.g_al,self.PRE))
            if len(g_al_input) != 0:
                self.g_al = float(g_al_input)
            g_be_input = input(' Lattice beta ({0:.{1}f})? '.format(self.g_be,self.PRE))
            if len(g_be_input) != 0:
                self.g_be = float(g_be_input)
            g_ga_input = input(' Lattice gam ({0:.{1}f})? '.format(self.g_ga,self.PRE))
            if len(g_ga_input) != 0:
                self.g_ga = float(g_ga_input)
        elif len(args) == 6:
            for value in args:
                if (type(value) in [int,float]) == False:
                    print ('\nInvalid argument: {0}\n'.format(value))
                    return
            self.g_aa, self.g_bb, self.g_cc, self.g_al, self.g_be, self.g_ga = args
        else:
            print ('\nUsage:  setlat()  or  setlat(a,b,c,alpha,beta,gam)\n')
            return
        g_sample_input = input('\nSample description: ({0})? '.format(self.g_sample))
        if len(g_sample_input) != 0:
            self.g_sample = g_sample_input
        print ('    -> Sample name set to {0}'.format(self.g_sample))
        # Calculate matrix U and B using new lattice parameters
        self.UB()
        # Get (H, K, L, ALPHA, BETA, AZIMUTH) using new lattice parameters
        self.wh_refresh()
    
    # Usage: setfrozen() or e.g., setfrozen(456)  or  setmode('045')
    def setfrozen(self, *args):
        """(), setfrozen(456) \t\t\t\t Set frozen of six-circle calculation. Choose which angles to freeze"""
        
        dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
        dic_pos = {0:self.F_TTH, 1:self.F_TH, 2:self.F_CHI, 3:self.F_PHI, 4:self.F_MU, 5:self.F_GAM, 6:self.F_OMEGA, 7:self.F_AZIMUTH, 8:self.F_ALPHA, 9:self.F_BETA}
        if len(args) == 0:
            print ('')
            print ('Current frozen: {0}'.format(self.g_frozen))
            g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(self.g_frozen)]
            print ('Current frozen angles:')
            print (' {0:>20}{1:>20}{2:>20}'.format(dic_ang[g_frozen_1],dic_ang[g_frozen_2],dic_ang[g_frozen_3]))
            print (' {0:>20.{3}f}{1:>20.{3}f}{2:>20.{3}f}'.format(dic_pos[g_frozen_1],dic_pos[g_frozen_2],dic_pos[g_frozen_3],self.PRE))
            print ('')
            print ('tth(0)  th(1)  chi(2)  phi(3)  mu(4)  gam(5)  omega(6)  azimuth(7)  alpha(8)  beta(9)')
            loop_flag = True
            while loop_flag == True:
                print ('')
                g_frozen_input = input('Select three frozen angles (a three-digit integer, e.g. 456): ')
                # If no input, use current setfrozen
                if len(g_frozen_input) == 0:
                    g_frozen_input = str(self.g_frozen)
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
        self.g_frozen = ''.join(sorted(list(g_frozen_input)))
        g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(self.g_frozen)]
        print ('')
        print ('Current frozen: {0}'.format(self.g_frozen))
        print ('Current frozen angles:')
        print (' {0:>20}{1:>20}{2:>20}'.format(dic_ang[g_frozen_1],dic_ang[g_frozen_2],dic_ang[g_frozen_3]))
        print (' {0:>20.{3}f}{1:>20.{3}f}{2:>20.{3}f}'.format(dic_pos[g_frozen_1],dic_pos[g_frozen_2],dic_pos[g_frozen_3],self.PRE))
        print ('Use freeze() command to change frozen values.')
        print ('')

    # Set positions of three frozen angles
    # Usage: freeze() or freeze(position1, position2, position3)
    def freeze(self, *args):
        """(), freeze(a1,a2,a3) \t\t\t\t\t Choose values for frozen angles (degrees)"""
        
        dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
        dic_pos = {0:self.F_TTH, 1:self.F_TH, 2:self.F_CHI, 3:self.F_PHI, 4:self.F_MU, 5:self.F_GAM, 6:self.F_OMEGA, 7:self.F_AZIMUTH, 8:self.F_ALPHA, 9:self.F_BETA}
        g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(self.g_frozen)]
        def set_freeze_positions(g_frozen_i,position_i):
            if g_frozen_i == 0:
                self.F_TTH = position_i
            elif g_frozen_i == 1:
                self.F_TH = position_i
            elif g_frozen_i == 2:
                self.F_CHI = position_i
            elif g_frozen_i == 3:
                self.F_PHI = position_i
            elif g_frozen_i == 4:
                self.F_MU = position_i
            elif g_frozen_i == 5:
                self.F_GAM = position_i
            elif g_frozen_i == 6:
                self.F_OMEGA = position_i
            elif g_frozen_i == 7:
                self.F_AZIMUTH = position_i
            elif g_frozen_i == 8:
                self.F_ALPHA = position_i
            elif g_frozen_i == 9:
                self.F_BETA = position_i
        if len(args) == 0:
            print('')
            position_1_input = input(' Freeze {0} ({1:.{2}f})? '.format(dic_ang[g_frozen_1], dic_pos[g_frozen_1], self.PRE))
            if len(position_1_input) != 0:
                position_1 = float(position_1_input)
                set_freeze_positions(g_frozen_1,position_1)
            position_2_input = input(' Freeze {0} ({1:.{2}f})? '.format(dic_ang[g_frozen_2], dic_pos[g_frozen_2], self.PRE))
            if len(position_2_input) != 0:
                position_2 = float(position_2_input)
                set_freeze_positions(g_frozen_2,position_2)
            position_3_input = input(' Freeze {0} ({1:.{2}f})? '.format(dic_ang[g_frozen_3], dic_pos[g_frozen_3], self.PRE))
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
        dic_pos = {0:self.F_TTH, 1:self.F_TH, 2:self.F_CHI, 3:self.F_PHI, 4:self.F_MU, 5:self.F_GAM, 6:self.F_OMEGA, 7:self.F_AZIMUTH, 8:self.F_ALPHA, 9:self.F_BETA}
        print ('')
        print ('Positions of frozen angles:')
        print ('{0:>10}{1:>10}{2:>10}'.format(dic_ang[g_frozen_1],dic_ang[g_frozen_2],dic_ang[g_frozen_3]))
        print ('{0:>10.{3}f}{1:>10.{3}f}{2:>10.{3}f}'.format(dic_pos[g_frozen_1],dic_pos[g_frozen_2],dic_pos[g_frozen_3],self.PRE))

    # Set azimuth reference vector
    # Usage: setaz() or setaz(H,K,L)
    def setaz(self, *args):
        """(), setaz(H,K,L) \t\t\t\t\t\t Sets surface normal/azimuthal reference"""
        
        if len(args) == 0:
            print ('\nEnter azimuthal reference H K L:')
            g_haz_input = input(' Azimuthal H ({0})? '.format(self.g_haz))
            if len(g_haz_input) != 0:
                self.g_haz = float(g_haz_input)
            g_kaz_input = input(' Azimuthal K ({0})? '.format(self.g_kaz))
            if len(g_kaz_input) != 0:
                self.g_kaz = float(g_kaz_input)        
            g_laz_input = input(' Azimuthal L ({0})? '.format(self.g_laz))
            if len(g_laz_input) != 0:
                self.g_laz = float(g_laz_input)
        elif len(args) == 3:
            for value in args:
                if (type(value) in [int,float]) == False:
                    print ('\nInvalid argument: {0}\n'.format(value))
                    return
            self.g_haz, self.g_kaz, self.g_laz = args
        else:
            print ('\nUsage:  setaz()  or  setaz(H,K,L)\n')
            return
        # Get (..., ALPHA, BETA, AZIMUTH) using new azimuth reference vector
        self.wh_refresh()

    # Positions of primary reflection are current positions
    # Usage: or0() or or0(H,K,L)
    def or0(self, *args):
        """(H,K,L), or1(H,K,L) \t\t\t\t\t\t Set H, K, L of primary/secondary reflection. Set primary,secondary (or0,or1) at present angle values"""

        self.g_lambda0 = self.LAMBDA
        self.g_u00, self.g_u01, self.g_u02, self.g_u03, self.g_u04, self.g_u05 = (self.TTH, self.TH, self.CHI, self.PHI, self.MU, self.GAM)
        if len(args) == 0:
            print ('\nEnter primary-reflection HKL coordinates:')
            g_h0_input = input(' H ({0:.{1}f})? '.format(self.g_h0,self.PRE))
            g_k0_input = input(' K ({0:.{1}f})? '.format(self.g_k0,self.PRE))
            g_l0_input = input(' L ({0:.{1}f})? '.format(self.g_l0,self.PRE))
            if len(g_h0_input) != 0:
                self.g_h0 = float(g_h0_input)
            if len(g_k0_input) != 0:
                self.g_k0 = float(g_k0_input)
            if len(g_l0_input) != 0:
                self.g_l0 = float(g_l0_input)
        elif len(args) == 3:
            for value in args:
                if (type(value) in [int,float]) == False:
                    print ('\nInvalid argument: {0}\n'.format(value))
                    return
            self.g_h0, self.g_k0, self.g_l0 = args
        else:
            print ('\nUsage:  or0()  or  or0(H,K,L)\n')
            return
        self.UB()
        self.wh_refresh()

    # Set H, K, L and positions of primary reflection
    def setor0(self):
        """() , setor1() \t\t\t\t\t\t Set primary,secondary (or0,or1) at entered angles"""

        self.g_lambda0 = self.LAMBDA
        print ('\nEnter primary-reflection angles:')
        g_u00_input = input(' Two Theta ({0:.{1}f})? '.format(self.g_u00,self.PRE))
        if len(g_u00_input) != 0:
            self.g_u00 = float(g_u00_input)
        g_u01_input = input(' Theta ({0:.{1}f})? '.format(self.g_u01,self.PRE))
        if len(g_u01_input) != 0:
            self.g_u01 = float(g_u01_input)
        g_u02_input = input(' Chi ({0:.{1}f})? '.format(self.g_u02,self.PRE))
        if len(g_u02_input) != 0:
            self.g_u02 = float(g_u02_input)
        g_u03_input = input(' Phi ({0:.{1}f})? '.format(self.g_u03,self.PRE))
        if len(g_u03_input) != 0:
            self.g_u03 = float(g_u03_input)
        g_u04_input = input(' Mu ({0:.{1}f})? '.format(self.g_u04,self.PRE))
        if len(g_u04_input) != 0:
            self.g_u04 = float(g_u04_input)
        g_u05_input = input(' Gam ({0:.{1}f})? '.format(self.g_u05,self.PRE))
        if len(g_u05_input) != 0:
            self.g_u05 = float(g_u05_input)
        print ('\nEnter primary-reflection HKL coordinates:')
        g_h0_input = input(' H ({0:.{1}f})? '.format(self.g_h0,self.PRE))
        if len(g_h0_input) != 0:
            self.g_h0 = float(g_h0_input)
        g_k0_input = input(' K ({0:.{1}f})? '.format(self.g_k0,self.PRE))
        if len(g_k0_input) != 0:
            self.g_k0 = float(g_k0_input)
        g_l0_input = input(' L ({0:.{1}f})? '.format(self.g_l0,self.PRE))
        if len(g_l0_input) != 0:
            self.g_l0 = float(g_l0_input)
        self.UB()
        self.wh_refresh()

    # Positions of secondary reflection are current positions
    # Usage: or1() or or1(H,K,L)
    def or1(self, *args):
        
        self.g_lambda1 = self.LAMBDA
        self.g_u10, self.g_u11, self.g_u12, self.g_u13, self.g_u14, self.g_u15 = (self.TTH, self.TH, self.CHI, self.PHI, self.MU, self.GAM)
        if len(args) == 0:
            print ('\nEnter secondary-reflection HKL coordinates:')
            g_h1_input = input(' H ({0:.{1}f})? '.format(self.g_h1,self.PRE))
            g_k1_input = input(' K ({0:.{1}f})? '.format(self.g_k1,self.PRE))
            g_l1_input = input(' L ({0:.{1}f})? '.format(self.g_l1,self.PRE))
            if len(g_h1_input) != 0:
                self.g_h1 = float(g_h1_input)
            if len(g_k1_input) != 0:
                self.g_k1 = float(g_k1_input)
            if len(g_l1_input) != 0:
                self.g_l1 = float(g_l1_input)
        elif len(args) == 3:
            for value in args:
                if (type(value) in [int,float]) == False:
                    print ('\nInvalid argument: {0}\n'.format(value))
                    return
            self.g_h1, self.g_k1, self.g_l1 = args
        else:
            print ('\nUsage:  or1()  or  or1(H,K,L)\n')
            return
        self.UB()
        self.wh_refresh()

    # Set H, K, L and positions of secondary reflection
    def setor1(self):
        
        self.g_lambda1 = self.LAMBDA
        print ('\nEnter secondary-reflection angles:')
        g_u10_input = input(' Two Theta ({0:.{1}f})? '.format(self.g_u10,self.PRE))
        if len(g_u10_input) != 0:
            self.g_u10 = float(g_u10_input)
        g_u11_input = input(' Theta ({0:.{1}f})? '.format(self.g_u11,self.PRE))
        if len(g_u11_input) != 0:
            self.g_u11 = float(g_u11_input)
        g_u12_input = input(' Chi ({0:.{1}f})? '.format(self.g_u12,self.PRE))
        if len(g_u12_input) != 0:
            self.g_u12 = float(g_u12_input)
        g_u13_input = input(' Phi ({0:.{1}f})? '.format(self.g_u13,self.PRE))
        if len(g_u13_input) != 0:
            self.g_u13 = float(g_u13_input)
        g_u14_input = input(' Mu ({0:.{1}f})? '.format(self.g_u14,self.PRE))
        if len(g_u14_input) != 0:
            self.g_u14 = float(g_u14_input)
        g_u15_input = input(' Gam ({0:.{1}f})? '.format(self.g_u15,self.PRE))
        if len(g_u15_input) != 0:
            self.g_u15 = float(g_u15_input)
        print ('\nEnter secondary-reflection HKL coordinates:')
        g_h1_input = input(' H ({0:.{1}f})? '.format(self.g_h1,self.PRE))
        if len(g_h1_input) != 0:
            self.g_h1 = float(g_h1_input)
        g_k1_input = input(' K ({0:.{1}f})? '.format(self.g_k1,self.PRE))
        if len(g_k1_input) != 0:
            self.g_k1 = float(g_k1_input)
        g_l1_input = input(' L ({0:.{1}f})? '.format(self.g_l1,self.PRE))
        if len(g_l1_input) != 0:
            self.g_l1 = float(g_l1_input)
        self.UB()
        self.wh_refresh()

    # Check consistency at or0 and or1
    def or_check(self):
        """() \t\t\t\t\t\t\t\t\t Check consistency of present or0 and or1 values"""
        
        # Record original conditions
        o_FLAG_WH = self.FLAG_WH
        o_tth, o_th, o_chi, o_phi, o_mu, o_gam = (self.TTH, self.TH, self.CHI, self.PHI, self.MU, self.GAM)
        self.FLAG_WH = False
        # Deviation of H,K,L at or0
        self.mv(tth=self.g_u00,th=self.g_u01,chi=self.g_u02,phi=self.g_u03,mu=self.g_u04,gam=self.g_u05)
        dH0 = self.H - self.g_h0
        dK0 = self.K - self.g_k0
        dL0 = self.L - self.g_l0
        # Deviation of H,K,L at or1
        self.mv(tth=self.g_u10,th=self.g_u11,chi=self.g_u12,phi=self.g_u13,mu=self.g_u14,gam=self.g_u15)
        dH1 = self.H - self.g_h1
        dK1 = self.K - self.g_k1
        dL1 = self.L - self.g_l1
        print ('Check self consistency {} '.format(self.g_sample),end='')
        print (' with ({0:.{3}f}, {1:.{3}f}, {2:.{3}f} '.format(self.g_aa,self.g_bb,self.g_cc,self.PRE),end='')
        print (' {1:.{3}f}, {2:.{3}f}, {2:.{3}f})'.format(self.g_al,self.g_be,self.g_ga,2))
        print ('   At or0 ({0:.{6}f}, {1:.{6}f}, {2:.{6}f}): dH = {3:.{6}f}, dK = {4:.{6}f}, dL = {5:.{6}f}'.format(self.g_h0,self.g_k0,self.g_l0,dH0,dK0,dL0,self.PRE))
        print ('   At or1 ({0:.{6}f}, {1:.{6}f}, {2:.{6}f}): dH = {3:.{6}f}, dK = {4:.{6}f}, dL = {5:.{6}f}'.format(self.g_h1,self.g_k1,self.g_l1,dH1,dK1,dL1,self.PRE))
        print ('')
        # Return to original conditions
        self.mv(tth=o_tth,th=o_th,chi=o_chi,phi=o_phi,mu=o_mu,gam=o_gam)
        self.FLAG_WH = o_FLAG_WH

    # Swap primary reflection and secondary reflection
    def or_swap(self):
        """() \t\t\t\t\t\t\t\t\t Swaps primary and secondary reference vectors/angles """

        or0 = (self.g_lambda0, self.g_u00, self.g_u01, self.g_u02, self.g_u03, self.g_u04, self.g_u05, self.g_h0, self.g_k0, self.g_l0)
        or1 = (self.g_lambda1, self.g_u10, self.g_u11, self.g_u12, self.g_u13, self.g_u14, self.g_u15, self.g_h1, self.g_k1, self.g_l1)
        self.g_lambda0, self.g_u00, self.g_u01, self.g_u02, self.g_u03, self.g_u04, self.g_u05, self.g_h0, self.g_k0, self.g_l0 = or1
        self.g_lambda1, self.g_u10, self.g_u11, self.g_u12, self.g_u13, self.g_u14, self.g_u15, self.g_h1, self.g_k1, self.g_l1 = or0
        self.UB()
        self.wh_refresh()

    # Redirecting to setfrozen
    def setmode(self):
        print ('Redirecting to setfrozen')
        self.setfrozen()

    # Print H, K, L, ALPHA, BETA, AZIMUTH, SA, Q
    def wh_refresh(self):

        self.OMEGA = self.TH - self.TTH/2
        # Scattering angle
        self.SA = 2*scbasic.thetaD_angle_2(self.TTH/2,self.MU,self.GAM)
        # Scattering vector
        self.ABSQ = scbasic.Q_length(self.LAMBDA,self.SA/2)
        # Azimuth reference vector
        az = np.array([[self.g_haz],[self.g_kaz],[self.g_laz]], dtype=float)
        hb = scbasic.CheckBr(self.LAMBDA,self.M_U,self.M_B,self.TTH/2,self.OMEGA,self.CHI,self.PHI,self.MU,self.GAM)
        self.H = hb[0][0]
        self.K = hb[1][0]
        self.L = hb[2][0]
        # When azimuth reference vector is surface normal, ALPHA and BETA are incident and exit angles
        self.AZIMUTH, self.ALPHA, self.BETA = scbasic.CheckPsiAlphainBetaout(self.M_U,self.M_B,az,self.TTH/2,self.OMEGA,self.CHI,self.PHI,self.MU,self.GAM)

    # Show current positions
    def wh(self):
        """()  (==wa()) \t\t\t\t\t\t\t\t Prints present HKL, etc.  (SA==Scattering Angle)"""
        print('')
        print ('H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(self.H,self.K,self.L,self.PRE))
        print ('|Q| = {0:.3f} nm-1  SA = {1:.{3}f} deg  at  LAMBDA = {2:.{4}f} A'.format(self.ABSQ*10,self.SA,self.LAMBDA,self.PRE,self.PRE+2))
        print ('AZ = ({0}, {1}, {2})  AZIMUTH = {3:.{6}f} deg  ALPHA = {4:.{6}f}  BETA = {5:.{6}f}'.format(self.g_haz,self.g_kaz,self.g_laz,self.AZIMUTH,self.ALPHA,self.BETA,self.PRE))
        print ('Omega = th-tth/2 = {0:.{1}f}'.format(self.OMEGA,self.PRE))
        print ('')
        strfmt = ('{:>'+str(self.PRE+6)+'}')*6
        strprt = ('tth','th','chi','phi','mu','gam')
        print (strfmt.format(*strprt))
        posfmt = ('{:>'+str(self.PRE+6)+'.'+str(self.PRE)+'f}')*6
        posprt = (self.TTH,self.TH,self.CHI,self.PHI,self.MU,self.GAM)
        print (posfmt.format(*posprt))
        print ('')

    # Redirecting to wh()
    def wa(self):
        print ('Redirecting to wh')
        self.wh()

    # Start printing positions after mv() or br()
    def wh_on(self):
        """(), wh_off() \t\t\t\t\t\t\t turns on/off wh() call at end of br() and mv() calls"""
        
        self.FLAG_WH = True
        print ('')
        print ('Will print positions after mv() or br()')
        print ('')

    # End printing positions after mv() or br()
    def wh_off(self):
        self.FLAG_WH = False
        print ('')
        print ('End printing positions after mv() or br()')
        print ('')

    # Change motors positions
    # Usage: mv(mu=?,gam=?,tth=?,th=?,chi=?,phi=?)
    def mv(self, **args):
        """(tth=?,th=?,chi=?,phi=?,mu=?,gam=?) \t\t Change angles"""
        
        if len(args) == 0:
            print ('\nUsage: mv(tth=?,th=?,chi=?,phi=?,mu=?,gam=?)\n')
            return

        for key in args.keys():
            if (key in self.mnemonics) == False:
                print ('Invalid motor mnemonic for mv: {0}\n'.format(key))
                print ('Valid motor mnemonic for mv: tth, th, chi, phi, mu, gam')
                return
        for value in args.values():
            if (type(value) in [int,float]) == False:
                print ('Invalid motor position for mv: {0}\n'.format(value))
                return
        if 'tth' in args.keys():
            self.TTH = args['tth']
        if 'th' in args.keys():
            self.TH = args['th']
        if 'chi' in args.keys():
            self.CHI = args['chi']
        if 'phi' in args.keys():
            self.PHI = args['phi']
        if 'mu' in args.keys():
            self.MU = args['mu']
        if 'gam' in args.keys():
            self.GAM = args['gam']
        # Get (H, K, L, SA, OMEGA, ALPHA, BETA, AZIMUTH, ABSQ) at new positions
        self.wh_refresh()
        if self.FLAG_WH == True:
            self.wh()

    # Set limits of positions
    # Usage:  setlm()  or  setlm(ltth=?,utth=?,lth=?,uth=?,lchi=?,uchi=?,lphi=?,uphi=?,lmu=?,umu=?,lgam=?,ugam=?,lalpha=?,ualpha=?,lbeta=?,ubeta=?)
    def setlm(self, **args):
        """() , showlm() \t\t\t\t\t\t\t Sets/Shows limits.  Useful selecting amoung multiple solutions for (H,K,L)"""

        if len(args) == 0:
            print ('\nSet limit of positions:')
            print ('')
            L_TTH_input = input(' Lower limit of tth ({0:.{1}f})? '.format(self.L_TTH,self.PRE))
            if len(L_TTH_input) != 0:
                self.L_TTH = float(L_TTH_input)
            U_TTH_input = input(' Upper limit of tth ({0:.{1}f})? '.format(self.U_TTH,self.PRE))
            if len(U_TTH_input) != 0:
                self.U_TTH = float(U_TTH_input)
            print ('')
            L_TH_input = input(' Lower limit of th ({0:.{1}f})? '.format(self.L_TH,self.PRE))
            if len(L_TH_input) != 0:
                self.L_TH = float(L_TH_input)
            U_TH_input = input(' Upper limit of th ({0:.{1}f})? '.format(self.U_TH,self.PRE))
            if len(U_TH_input) != 0:
                self.U_TH = float(U_TH_input)
            print ('')
            L_CHI_input = input(' Lower limit of chi ({0:.{1}f})? '.format(self.L_CHI,self.PRE))
            if len(L_CHI_input) != 0:
                self.L_CHI = float(L_CHI_input)
            U_CHI_input = input(' Upper limit of chi ({0:.{1}f})? '.format(self.U_CHI,self.PRE))
            if len(U_CHI_input) != 0:
                self.U_CHI = float(U_CHI_input)
            print ('')
            L_PHI_input = input(' Lower limit of phi ({0:.{1}f})? '.format(self.L_PHI,self.PRE))
            if len(L_PHI_input) != 0:
                self.L_PHI = float(L_PHI_input)
            U_PHI_input = input(' Upper limit of phi ({0:.{1}f})? '.format(self.U_PHI,self.PRE))
            if len(U_PHI_input) != 0:
                self.U_PHI = float(U_PHI_input)
            print ('')
            L_MU_input = input(' Lower limit of mu ({0:.{1}f})? '.format(self.L_MU,self.PRE))
            if len(L_MU_input) != 0:
                self.L_MU = float(L_MU_input)
            U_MU_input = input(' Upper limit of mu ({0:.{1}f})? '.format(self.U_MU,self.PRE))
            if len(U_MU_input) != 0:
                self.U_MU = float(U_MU_input)
            print ('')
            L_GAM_input = input(' Lower limit of gam ({0:.{1}f})? '.format(self.L_GAM,self.PRE))
            if len(L_GAM_input) != 0:
                self.L_GAM = float(L_GAM_input)
            U_GAM_input = input(' Upper limit of gam ({0:.{1}f})? '.format(self.U_GAM,self.PRE))
            if len(U_GAM_input) != 0:
                self.U_GAM = float(U_GAM_input)
            print ('')
            L_ALPHA_input = input(' Lower limit of alpha ({0:.{1}f})? '.format(self.L_ALPHA,self.PRE))
            if len(L_ALPHA_input) != 0:
                self.L_ALPHA = float(L_ALPHA_input)
            U_ALPHA_input = input(' Upper limit of alpha ({0:.{1}f})? '.format(self.U_ALPHA,self.PRE))
            if len(U_ALPHA_input) != 0:
                self.U_ALPHA = float(U_ALPHA_input)
            print ('')
            L_BETA_input = input(' Lower limit of beta ({0:.{1}f})? '.format(self.L_BETA,self.PRE))
            if len(L_BETA_input) != 0:
                self.L_BETA = float(L_BETA_input)
            U_BETA_input = input(' Upper limit of beta ({0:.{1}f})? '.format(self.U_BETA,self.PRE))
            if len(U_BETA_input) != 0:
                self.U_BETA = float(U_BETA_input)
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
                self.L_TTH = args['ltth']
            if 'utth' in args.keys():
                self.U_TTH = args['utth']
            if 'lth' in args.keys():
                self.L_TH = args['lth']
            if 'uth' in args.keys():
                self.U_TH = args['uth']
            if 'lchi' in args.keys():
                self.L_CHI = args['lchi']
            if 'uchi' in args.keys():
                self.U_CHI = args['uchi']
            if 'lphi' in args.keys():
                self.L_PHI = args['lphi']
            if 'uphi' in args.keys():
                self.U_PHI = args['uphi']
            if 'lmu' in args.keys():
                self.L_MU = args['lmu']
            if 'umu' in args.keys():
                self.U_MU = args['umu']
            if 'lgam' in args.keys():
                self.L_GAM = args['lgam']
            if 'ugam' in args.keys():
                self.U_GAM = args['ugam']
            if 'lalpha' in args.keys():
                self.L_ALPHA = args['lalpha']
            if 'lalpha' in args.keys():
                self.U_ALPHA = args['ualpha']
            if 'lbeta' in args.keys():
                self.L_BETA = args['lbeta']
            if 'ubeta' in args.keys():
                self.U_BETA = args['ubeta']
            self.showmainlm()

    def showlm(self):
        print("\n  Motor limits are now",end='')
        print("    tth: {} {} ".format(self.L_TTH,self.U_TTH), end='')
        print("     th: {} {} ".format(self.L_TH,self.U_TH), end='')
        print("    chi: {} {} ".format(self.L_CHI,self.U_CHI), end='')
        print("    phi: {} {} ".format(self.L_PHI,self.U_PHI), end='')
        print("     mu: {} {} ".format(self.L_MU,self.U_MU), end='')
        print("    gam: {} {} ".format(self.L_GAM,self.U_GAM), end='')
        print("  alpha: {} {} ".format(self.L_ALPHA,self.U_ALPHA), end='')
        print("   beta: {} {} ".format(self.L_BETA,self.U_BETA), end='')
        print("\n")

    def showmainlm(self):
        print("\n  Main motor limits are now",end='')
        print("    tth: {} {} ".format(self.L_TTH,self.U_TTH), end='')
        print("     th: {} {} ".format(self.L_TH,self.U_TH), end='')
        print("    chi: {} {} ".format(self.L_CHI,self.U_CHI), end='')
        print("    phi: {} {} ".format(self.L_PHI,self.U_PHI), end='')
        print("\n")
    
    # Set all limits as -180 to 180
    def setlm_clear(self):
        """() \t\t\t\t\t\t\t\t Sets default (+-180 deg) limits on all angles"""
        
        self.L_TTH = -180.0; self.U_TTH = 180.0
        self.L_TH = -180.0; self.U_TH = 180.0
        self.L_CHI = -180.0; self.U_CHI = 180.0 
        self.L_PHI = -180.0; self.U_PHI = 180.0
        self.L_MU = -180.0; self.U_MU = 180.0
        self.L_GAM = -180.0; self.U_GAM = 180.0
        self.L_ALPHA = -180.0; self.U_ALPHA = 180.0 
        self.L_BETA = -180.0; self.U_BETA = 180.0
        print ('')
        print ('Now all limits are set as -180 to 180')
        print ('')

    # Calculate positions in specified frozen
    # Show the 1st set of positions in preset limits
    # Usage:  ca(H,K,L)
    def ca(self, *args):
        """(H,K,L) \t\t\t\t\t\t\t\t\t Calculate angles for a given (H,K,L) """
        
        if len(args) != 3:
            print ('Usage:  ca(H,K,L)')
            return
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        caH, caK, caL = args
        flag, pos = self.ca_s(caH, caK, caL)
        if flag == False:
            dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
            dic_pos = {0:self.F_TTH, 1:self.F_TH, 2:self.F_CHI, 3:self.F_PHI, 4:self.F_MU, 5:self.F_GAM, 6:self.F_OMEGA, 7:self.F_AZIMUTH, 8:self.F_ALPHA, 9:self.F_BETA}
            g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(self.g_frozen)]
            print ('Error: Impossible reflection within current limits for frozen: ', end='')
            freprt = (self.g_frozen,dic_ang[g_frozen_1],dic_pos[g_frozen_1],dic_ang[g_frozen_2],dic_pos[g_frozen_2],dic_ang[g_frozen_3],dic_pos[g_frozen_3])
            print ('{1}={2:.{7}f} {3}={4:.{7}f} {5}={6:.{7}f}'.format(*freprt,self.PRE))
            return
        else:
            caTTH, caTH, caCHI, caPHI, caMU, caGAM, caSA, caOMEGA, caAZIMUTH, caALPHA, caBETA = pos[0]
            caABSQ = scbasic.Q_length(self.LAMBDA,abs(caSA/2))
            print('')
            print ('H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(caH,caK,caL,self.PRE))
            print ('|Q| = {0:.3f} nm-1  SA = {1:.{3}f} deg  at  LAMBDA = {2:.{4}f} A'.format(caABSQ*10,caSA,self.LAMBDA,self.PRE,self.PRE+2))
            print ('AZ = ({0}, {1}, {2})  AZIMUTH = {3:.{6}f} deg  ALPHA = {4:.{6}f}  BETA = {5:.{6}f}'.format(self.g_haz,self.g_kaz,self.g_laz,caAZIMUTH,caALPHA,caBETA,self.PRE))
            print ('Omega = th-tth/2 = {0:.{1}f}'.format(caOMEGA,self.PRE))
            print ('')
            strfmt = ('{:>'+str(self.PRE+6)+'}')*6
            strprt = ('tth','th','chi','phi','mu','gam')
            print (strfmt.format(*strprt))
            posfmt = ('{:>'+str(self.PRE+6)+'.'+str(self.PRE)+'f}')*6
            posprt = (caTTH,caTH,caCHI,caPHI,caMU,caGAM)
            print (posfmt.format(*posprt))
            print ('')
            print ('Command (sixcircle):  ', end='')
            print ('mv (tth={0:.{6}f}, th={1:.{6}f}, chi={2:.{6}f}, phi={3:.{6}f}, mu={4:.{6}f}, gam={5:.{6}f})'.format(*posprt,self.PRE))
            print ('Command (BL43LXU):    ', end='')
            print ('mv tth {0:.{6}f} th {1:.{6}f} chi {2:.{6}f} phi {3:.{6}f}'.format(*posprt,self.PRE))
            print ('')
    
    # Move to H, K, L    Using first set of positions from ca_s()
    # Usage: br(H,K,L)
    def br(self, *args):
        """(H,K,L) \t\t\t\t\t\t\t\t\t Move to given Q=(H,K,L)"""
        
        if len(args) != 3:
            print ('Usage:  br(H,K,L)')
            return
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        goH, goK, goL = args
        flag, pos = self.ca_s(goH, goK, goL)
        dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
        dic_pos = {0:self.F_TTH, 1:self.F_TH, 2:self.F_CHI, 3:self.F_PHI, 4:self.F_MU, 5:self.F_GAM, 6:self.F_OMEGA, 7:self.F_AZIMUTH, 8:self.F_ALPHA, 9:self.F_BETA}
        g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(self.g_frozen)]
        if flag == False:
            print ('Error: Impossible reflection within current limits for frozen {0}: ', end='')
            freprt = (self.g_frozen,dic_ang[g_frozen_1],dic_pos[g_frozen_1],dic_ang[g_frozen_2],dic_pos[g_frozen_2],dic_ang[g_frozen_3],dic_pos[g_frozen_3])
            print ('freeze {1}={2:.{7}f} {3}={4:.{7}f} {5}={6:.{7}f}'.format(*freprt,self.PRE))
            return
        goprt = (goH,goK,goL,dic_ang[g_frozen_1],dic_ang[g_frozen_2],dic_ang[g_frozen_3],dic_pos[g_frozen_1],dic_pos[g_frozen_2],dic_pos[g_frozen_3])
        print ('Moving to ({0}, {1}, {2}) with ({3}, {4}, {5}) frozen at ({6:.{9}f}, {7:.{9}f}, {8:.{9}f})'.format(*goprt,self.PRE))
        go_tth,go_th,go_chi,go_phi,go_mu,go_gam,go_sa,go_omega,go_azimuth,go_alpha,go_beta = pos[0]
        self.mv(tth=go_tth,th=go_th,chi=go_chi,phi=go_phi,mu=go_mu,gam=go_gam)
    
    # Calculate positions in specified frozen
    # _s: silent, calculation in background
    # return to: flag, pos
    # flag: True (False) when calculation is successful (failed: no solution)
    # pos: N set of positions, [[tth1, th1, chi1, phi1, mu1, gam1, sa1, omega1, azimuth1, alpha1, beta1],[tth2,...],...,[tthN,...]]
    def ca_s(self, caH, caK, caL):
        
        # Azimuth reference vector
        az = np.array([[self.g_haz],[self.g_kaz],[self.g_laz]], dtype=float)
        # Target H, K, L
        ca_hb = np.array([[caH],[caK],[caL]], dtype=float)
        # Pre-specified frozen and positions of frozen angles
        dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
        dic_pos = {0:self.F_TTH, 1:self.F_TH, 2:self.F_CHI, 3:self.F_PHI, 4:self.F_MU, 5:self.F_GAM, 6:self.F_OMEGA, 7:self.F_AZIMUTH, 8:self.F_ALPHA, 9:self.F_BETA}
        g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(self.g_frozen)]
        # Prepare for calculation. Put in positions of specified frozen angles. For unfrozen angles to be calculated: 'x'
        dic = {'tth':'x','th':'x','chi':'x','phi':'x','mu':'x','gam':'x','omega':'x','azimuth':'x','alpha':'x','beta':'x','theta':'x'}
        dic[dic_ang[g_frozen_1]] = dic_pos[g_frozen_1]
        dic[dic_ang[g_frozen_2]] = dic_pos[g_frozen_2]
        dic[dic_ang[g_frozen_3]] = dic_pos[g_frozen_3]
        # If freeze tth, use theta = tth/2 in calculation. Note theta != th
        g_frozen_str = str(self.g_frozen)
        if g_frozen_str.count('0') == 1:
            dic['theta'] = dic['tth']/2
        flag = False
        pos = []
        # In general: neither azimuth, alpha, beta is frozen
        if g_frozen_str.count('7') + g_frozen_str.count('8') + g_frozen_str.count('9') == 0:
            flagp, angles = scbasic.SC_angles(self.LAMBDA,self.M_U,self.M_B,ca_hb,dic['th'],dic['theta'],dic['omega'],dic['chi'],dic['phi'],dic['mu'],dic['gam'])
        else:
            # Azimuth frozen
            if g_frozen_str.count('7') == 1:
                flagp, angles = scbasic.SC_angles_fix_psi(self.LAMBDA,self.M_U,self.M_B,ca_hb,az,dic['azimuth'],dic['th'],dic['theta'],dic['omega'],dic['chi'],dic['phi'],dic['mu'],dic['gam'])
            # Alpha frozen or Beta frozen
            else:
                flagp, angles = scbasic.SC_angles_fix_ab(self.LAMBDA,self.M_U,self.M_B,ca_hb,az,dic['alpha'],dic['beta'],dic['th'],dic['theta'],dic['omega'],dic['chi'],dic['phi'],dic['mu'],dic['gam'])
        if flagp == True:
            # Current angles: [[theta,omega,chi,phi,mu,gam],...]
            # Sort out angles to [[tth,th,chi,phi,mu,gam],...]
            N, angles = scbasic.motors(angles)
            for element in angles:
                caTTH, caTH, caCHI, caPHI, caMU, caGAM = element
                caSA = 2*scbasic.thetaD_angle_2(caTTH/2,caMU,caGAM)
                caOMEGA = caTH - caTTH/2
                # Calculate corresponding azimuth, alpha, beta
                caAZIMUTH, caALPHA, caBETA = scbasic.CheckPsiAlphainBetaout(self.M_U,self.M_B,az,caTTH/2,caOMEGA,caCHI,caPHI,caMU,caGAM)
                if self.L_TTH <= caTTH < self.U_TTH and self.L_TH <= caTH < self.U_TH and self.L_CHI <= caCHI < self.U_CHI and self.L_PHI <= caPHI < self.U_PHI:
                    if self.L_MU <= caMU < self.U_MU and self.L_GAM <= caGAM < self.U_GAM:
                        if self.L_ALPHA <= caALPHA < self.U_ALPHA and self.L_BETA <= caBETA < self.U_BETA:
                            pos.append([caTTH,caTH,caCHI,caPHI,caMU,caGAM,caSA,caOMEGA,caAZIMUTH,caALPHA,caBETA])
        if len(pos) > 0:
            flag = True
        return flag, pos
    
    # Calculate positions in specified frozen
    # _a: all, show all sets of positions in preset limits
    # Usage:  ca_a(H,K,L)
    def ca_a(self, *args):
        """(H,K,L), ca_s(H,K,L) \t\t\t\t\t Calculates (_s=silently) angles for (H,K,L) within limits from setlm"""
        
        if len(args) != 3:
            print ('Usage:  ca_a(H,K,L)')
            return
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        caH, caK, caL = args
        flag, pos = self.ca_s(caH, caK, caL)
        if flag == False:
            dic_ang = {0:'tth', 1:'th', 2:'chi', 3:'phi', 4:'mu', 5:'gam', 6:'omega', 7:'azimuth', 8:'alpha', 9:'beta'}
            dic_pos = {0:self.F_TTH, 1:self.F_TH, 2:self.F_CHI, 3:self.F_PHI, 4:self.F_MU, 5:self.F_GAM, 6:self.F_OMEGA, 7:self.F_AZIMUTH, 8:self.F_ALPHA, 9:self.F_BETA}
            g_frozen_1, g_frozen_2, g_frozen_3 = [int(i) for i in list(self.g_frozen)]
            print ('Error: Impossible reflection within current limits for frozen {0}: ', end='')
            freprt = (self.g_frozen,dic_ang[g_frozen_1],dic_pos[g_frozen_1],dic_ang[g_frozen_2],dic_pos[g_frozen_2],dic_ang[g_frozen_3],dic_pos[g_frozen_3])
            print ('freeze {1}={2:.{7}f} {3}={4:.{7}f} {5}={6:.{7}f}'.format(*freprt,self.PRE))
            return
        else:
            caSA = pos[0][6]
            caABSQ = scbasic.Q_length(self.LAMBDA,abs(caSA/2))
            print ('')
            print ('Calculated Positions:')
            print ('')
            print ('H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(caH,caK,caL,self.PRE))
            print ('|Q| = {0:.3f} nm-1  SA = {1:.{3}f} deg  at  LAMBDA = {2:.{4}f} A'.format(caABSQ*10,caSA,self.LAMBDA,self.PRE,self.PRE+2))
            print ('')
            strfmt = '    '+('{:>'+str(self.PRE+6)+'}')*11
            strprt = ('tth','th','chi','phi','mu','gam','sa','omega','azimuth','alpha','beta')
            print (strfmt.format(*strprt))
            posfmt = '{:>4}'+('{:>'+str(self.PRE+6)+'.'+str(self.PRE)+'f}')*11
            for i in range(0,len(pos)):
                print (posfmt.format(i,*pos[i],self.PRE))
            print ('')
            print ('Command (sixcircle):')
            for i in range(0,len(pos)):
                posprt = [pos[i][0],pos[i][1],pos[i][2],pos[i][3],pos[i][4],pos[i][5]]
                print ('{0:>4}  mv (tth={1:.{7}f}, th={2:.{7}f}, chi={3:.{7}f}, phi={4:.{7}f}, mu={5:.{7}f}, gam={6:.{7}f})'.format(i,*posprt,self.PRE))        
            print ('')

    # Check where the limits of alpha and beta are for H, K, L
    # _s: silent, calculation in background
    # Usage:  flag, min, max = wmab_s(H, K, L)
    def wmab_s(self, caH, caK, caL):
        # Azimuth reference vector
        az = np.array([[self.g_haz],[self.g_kaz],[self.g_laz]], dtype=float)
        # Target H, K, L
        ca_hb = np.array([[caH],[caK],[caL]], dtype=float)
        flag, minab, maxab = scbasic.CheckRangeAlphainBetaout(self.LAMBDA,self.M_B,ca_hb,az)
        return flag, minab, maxab

    # Check where the limits of alpha and beta are for H, K, L
    # Usage:  wmab(H, K, L)
    def wmab(self, *args):
        """(H,K,L), wmab_s(H,K,L) \t\t\t\t\t Check the limits on alpha and beta for given (H,K,L)"""

        if len(args) != 3:
            print ('Usage:  wmab(H,K,L)')
            return
        for value in args:
            if (type(value) in [int,float]) == False:
                print ('\nInvalid argument: {0}\n'.format(value))
                return
        caH, caK, caL = args
        flag, minab, maxab = self.wmab_s(caH, caK, caL)
        if flag == False:
            print ('')
            print ('Error: Impossible reflection.')
            return
        print ('')
        print ('Limits of ALPHA and BETA for {0}, {1}, {2}'.format(caH, caK, caL))
        print ('')
        print ('{0:>10}{1:>10}'.format('Min','Max'))
        print ('{0:>10.{2}f}{1:>10.{2}f}'.format(minab,maxab,self.PRE))
        print ('')

    # Set output precision
    # Usage:  setprecision()  or setprecision(n)
    # Default n = 4
    def setprecision(self, *args):
        """(n) \t\t\t\t\t\t\t Sets decimal precision for output (default n=4)"""

        if len(args) == 0:
            PRE_input =  input('\nOutput precision ({0})? '.format(self.PRE))
            if len(PRE_input) != 0:
                self.PRE = int(PRE_input)
        elif len(args) == 1:
            if type(args[0]) != int:
                print ('\nInvalid argument: {0}\n'.format(args[0]))
                return
            if args[0] < 0:
                print ('\nInvalid argument: {0}\n'.format(args[0]))
                return          
            self.PRE = int(args[0])
        else:
            print ('\nUsage:  setprecision()  or  setprecision(n)\n')
            return
        print ('')
        print ('Output precision set to {0}'.format(self.PRE))
        print ('')
        
    def help(self):
        self.__doc__()