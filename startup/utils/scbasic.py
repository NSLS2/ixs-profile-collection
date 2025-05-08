# scbasic: Modules for six-circle calculation
# Required by sixcircle

# Developed with Python 3.7.3, numpy 1.16.4, scipy 1.3.0
# Materials Dynamics Laboratory, RIKEN SPring-8 Center
# Ver 1.5, December 2020

# Authors: Wenyang Zhao & Alfred Q. R. Baron
# Contact: baron@spring8.or.jp

# If this code is used independently of work at BL43, please reference the writeup ("Open-source Python software for six-circle diffraction with an inelastic x-ray scattering (IXS) spectrometer" by Wenyang ZHAO and Alfred Q.R. BARON, unpublished, available at https://beamline.harima.riken.jp/bl43lxu)

# In typical use related to experimental work at BL43LXU, it is enough to reference a BL publication.  However, if this code is used extensively, then please include a separate reference as above.


import math
import numpy as np
import random
from scipy.optimize import least_squares

print("Importing scbasic.py containing definitions for basic diffraction calculations")
print("   Includes: XX_matrix, theta_angle, CheckBr, SC_angles, etc")
print("")
PI = np.pi
# Matrix B
# vc (vector in cartesian space) = B dot vb (vector HKL in reciprocal space)
def B_matrix(a1,a2,a3,alpha1d,alpha2d,alpha3d):
    # a1, a2, a3, alpha1, alpha2, alpha3: lattice parameters in direct space
    # b1, b2, b3, beta1, beta2, beta3: lattice parameters in reciprocal space
    # Note: *d: unit degree; *r: unit rad
    if sum([alpha1d,alpha2d,alpha3d]) <= 2*max(alpha1d,alpha2d,alpha3d):
        # Invalid lattice parameters in direct space
        flag = False
        return flag, 'non', 'non', 'non', 'non', 'non', 'non', 'non'
    alpha1r = math.radians(alpha1d)
    alpha2r = math.radians(alpha2d)
    alpha3r = math.radians(alpha3d)
    cosalpha1 = math.cos(alpha1r)
    cosalpha2 = math.cos(alpha2r)
    cosalpha3 = math.cos(alpha3r)
    sinalpha1 = math.sin(alpha1r)
    sinalpha2 = math.sin(alpha2r)
    sinalpha3 = math.sin(alpha3r)
    # Volume of unit cell
    va = a1 * a2 * a3 * (1 - cosalpha1**2 - cosalpha2**2 - cosalpha3**2 + 2*cosalpha1*cosalpha2*cosalpha3) ** 0.5
    b1 =  2 * PI * a2 * a3 * sinalpha1 / va
    b2 =  2 * PI * a1 * a3 * sinalpha2 / va
    b3 =  2 * PI * a1 * a2 * sinalpha3 / va
    cosbeta1 = (cosalpha2 * cosalpha3 - cosalpha1) / (sinalpha2 * sinalpha3)
    cosbeta2 = (cosalpha1 * cosalpha3 - cosalpha2) / (sinalpha1 * sinalpha3)
    cosbeta3 = (cosalpha1 * cosalpha2 - cosalpha3) / (sinalpha1 * sinalpha2)
    sinbeta1 = (1 - cosbeta1**2) ** 0.5
    sinbeta2 = (1 - cosbeta2**2) ** 0.5
    sinbeta3 = (1 - cosbeta3**2) ** 0.5
    beta1d = math.degrees(math.acos(cosbeta1))
    beta2d = math.degrees(math.acos(cosbeta2))
    beta3d = math.degrees(math.acos(cosbeta3))
    B11 = b1
    B12 = b2 * cosbeta3
    B13 = b3 * cosbeta2
    B21 = 0
    B22 = b2 * sinbeta3
    B23 = - b3 * sinbeta2 * cosalpha1
    B31 = 0
    B32 = 0
    B33 = 2 * PI / a3
    B = np.array([[B11,B12,B13], [B21,B22,B23], [B31,B32,B33]], dtype=float)
    flag = True
    return flag, B, b1, b2, b3, beta1d, beta2d, beta3d

# Matrix PHI
# vchi (vector in chi space) = PHI dot vphi (vector in phi space)
def PHI_matrix(phid):
    phir = math.radians(phid)
    cosphi = math.cos(phir)
    sinphi = math.sin(phir)
    PHI11 = cosphi
    PHI12 = - sinphi
    PHI13 = 0
    PHI21 = sinphi
    PHI22 = cosphi
    PHI23 = 0
    PHI31 = 0
    PHI32 = 0
    PHI33 = 1
    PHI = np.array([[PHI11,PHI12,PHI13], [PHI21,PHI22,PHI23], [PHI31,PHI32,PHI33]], dtype=float)
    return PHI

# Matrix CHI
# vomega (vector in omega space) = CHI dot vchi (vector in chi space)
def CHI_matrix(chid):
    chir = math.radians(chid)
    coschi = math.cos(chir)
    sinchi = math.sin(chir)
    CHI11 = coschi
    CHI12 = 0
    CHI13 = - sinchi
    CHI21 = 0
    CHI22 = 1
    CHI23 = 0
    CHI31 = sinchi
    CHI32 = 0
    CHI33 = coschi
    CHI = np.array([[CHI11,CHI12,CHI13], [CHI21,CHI22,CHI23], [CHI31,CHI32,CHI33]], dtype=float)
    return CHI

# Matrix OMEGA
# vtheta (vector in theta space) = OMEGA dot vomega (vector in omega space)
def OMEGA_matrix(omegad):
    omegar = math.radians(omegad)
    cosomega = math.cos(omegar)
    sinomega = math.sin(omegar)
    OMEGA11 = cosomega
    OMEGA12 = - sinomega
    OMEGA13 = 0
    OMEGA21 = sinomega
    OMEGA22 = cosomega
    OMEGA23 = 0
    OMEGA31 = 0
    OMEGA32 = 0
    OMEGA33 = 1
    OMEGA = np.array([[OMEGA11,OMEGA12,OMEGA13], [OMEGA21,OMEGA22,OMEGA23], [OMEGA31,OMEGA32,OMEGA33]], dtype=float)
    return OMEGA

# Matrix THETA
# vl (vector in laboratory space) = THETA dot vtheta (vector in theta space)
# theta = tth/2
# thetaD = deflection_angle/2
# omega = th- tth/2
def THETA_matrix(thetad):
    thetar = math.radians(thetad)
    costheta = math.cos(thetar)
    sintheta = math.sin(thetar)
    THETA11 = costheta
    THETA12 = - sintheta
    THETA13 = 0
    THETA21 = sintheta
    THETA22 = costheta
    THETA23 = 0
    THETA31 = 0
    THETA32 = 0
    THETA33 = 1
    THETA = np.array([[THETA11,THETA12,THETA13], [THETA21,THETA22,THETA23], [THETA31,THETA32,THETA33]], dtype=float)
    return THETA

# Matrix PSI
# psi: azimuth angle
# vd (vector in diffraction-plane space) = PSI dot vpsi (vector in psi space)
def PSI_matrix(psid):
    psir = math.radians(psid)
    cospsi = math.cos(psir)
    sinpsi = math.sin(psir)
    PSI11 = 1
    PSI12 = 0
    PSI13 = 0
    PSI21 = 0
    PSI22 = cospsi
    PSI23 = - sinpsi
    PSI31 = 0
    PSI32 = sinpsi
    PSI33 = cospsi
    PSI = np.array([[PSI11, PSI12, PSI13], [PSI21, PSI22, PSI23] ,[PSI31, PSI32, PSI33]], dtype=float)
    return PSI

# Matrix D
# vl (vector in laboratory space) = D dot vd (vector in diffraction-plane space)
def D_matrix(thetad,mud,gammad):
    twothetad = thetad * 2
    twothetar = math.radians(twothetad)
    thetar = math.radians(thetad)
    mur = math.radians(mud)
    gammar = math.radians(gammad)
    sintwotheta = math.sin(twothetar)
    costwotheta = math.cos(twothetar)
    cosmu = math.cos(mur)
    sinmu = math.sin(mur)
    singamma = math.sin(gammar)
    cosgamma = math.cos(gammar)
    upl1 = 0
    upl2 = - cosmu
    upl3 = - sinmu
    upl = np.array([[upl1], [upl2], [upl3]], dtype=float)
    udl1 = sintwotheta * cosgamma
    udl2 = - costwotheta * cosgamma
    udl3 = - singamma
    udl = np.array([[udl1], [udl2], [udl3]], dtype=float)
    sl1 = udl1 - upl1
    sl2 = udl2 - upl2
    sl3 = udl3 - upl3
    sl = np.array([[sl1], [sl2], [sl3]], dtype=float)
    u1l = sl
    u3l = np.cross(u1l.flatten(), ((-1)*upl).flatten()).reshape(3,1)
    u2l = np.cross(u3l.flatten(), u1l.flatten()).reshape(3,1)
    if np.linalg.norm(sl) == 0:
        D = np.array([[1,0,0],[0,1,0],[0,0,1]], dtype=float)
    else:
        u1l = u1l/np.linalg.norm(u1l)
        u2l = u2l/np.linalg.norm(u2l)
        u3l = u3l/np.linalg.norm(u3l)
        D = np.hstack((u1l, u2l, u3l))
    return D

# uphi: unit scattering vector in phi space
def uphi_vector(thetad,omegad,chid,phid,mud,gammad):
    THETA_I = np.transpose(THETA_matrix(thetad))
    OMEGA_I = np.transpose(OMEGA_matrix(omegad))
    CHI_I = np.transpose(CHI_matrix(chid))
    PHI_I = np.transpose(PHI_matrix(phid))
    D = D_matrix(thetad,mud,gammad)
    ul = np.array([[1], [0], [0]], dtype=float)
    uphi = np.dot(np.dot(np.dot(np.dot(np.dot(PHI_I,CHI_I),OMEGA_I),THETA_I),D),ul)
    return uphi

# uc: unit vector of HKL in cartesian space
# hb = np.array([[H],[K],[L]], dtype=float)
def uc_vector(B,hb):
    uc = np.dot(B,hb)
    uc = uc / np.linalg.norm(uc)
    return uc

# Matrix U: orientation matrx
# vphi (vector in phi space) = U dot vc (vector in cartesian space)
# u1phi, u1c: primary reflection
# u2phi, u2c: secondary reflection
def U_matrix(u1phi, u2phi, u1c, u2c):
    U = np.array([[0,0,0],[0,0,0],[0,0,0]], dtype=float)
    errcode = 0
    t1c = u1c
    t3c = np.cross(t1c.flatten(), u2c.flatten()).reshape(3,1)
    t2c = np.cross(t3c.flatten(), t1c.flatten()).reshape(3,1)
    if abs(np.linalg.norm(t1c))<=1e-10 or abs(np.linalg.norm(t2c))<1e-10 or abs(np.linalg.norm(t3c))<1e-10:
        # Can't find orientation matrix:  Reflections are parallel.
        errcode = 1
        return errcode, U
    else:
        t1c = t1c / np.linalg.norm(t1c)
        t2c = t2c / np.linalg.norm(t2c)
        t3c = t3c / np.linalg.norm(t3c)
    t1phi = u1phi
    t3phi = np.cross(t1phi.flatten(), u2phi.flatten()).reshape(3,1)
    t2phi = np.cross(t3phi.flatten(), t1phi.flatten()).reshape(3,1)
    if abs(np.linalg.norm(t1phi))<1e-10 or abs(np.linalg.norm(t2phi))<1e-10 or abs(np.linalg.norm(t3phi))<1e-10:
        # Can't find orientation matrix:  Reflections (by angles) are parallel.
        errcode = 2
        return errcode, U
    else:
        t1phi = t1phi / np.linalg.norm(t1phi)
        t2phi = t2phi / np.linalg.norm(t2phi)
        t3phi = t3phi / np.linalg.norm(t3phi)
    Tc = np.hstack((t1c,t2c,t3c))
    Tphi = np.hstack((t1phi,t2phi,t3phi))
    U = np.dot(Tphi, np.transpose(Tc))
    return errcode, U

# Calculate half-deflection-angle (thetaD) using Bragg's law
def thetaD_angle(xlambda,B,hb):
    flag = True
    thetaDd = 0
    hc = np.dot(B,hb)
    sinthetaD = np.linalg.norm(hc)*xlambda/(4*PI)
    if abs(sinthetaD) > 1:
        # Impossible reflection HKL with present wavelength
        # return thetaD zero
        flag = False
        return flag, thetaDd
    else:
        thetaDr = math.asin(sinthetaD)
        thetaDd = math.degrees(thetaDr)
    return flag, thetaDd

# Calculate half-deflection-angle (thetaD) using theta, mu, gamma
def thetaD_angle_2(thetad,mud,gammad):
    twothetar = math.radians(thetad*2)
    costwotheta = math.cos(twothetar)
    mur = math.radians(mud)
    gammar = math.radians(gammad)
    sinmu = math.sin(mur)
    cosmu = math.cos(mur)
    singamma = math.sin(gammar)
    cosgamma = math.cos(gammar)
    sinthetaD = (2 - 2*sinmu*singamma - 2*cosmu*cosgamma*costwotheta) ** 0.5 / 2
    thetaDr = math.asin(sinthetaD)
    thetaDd = math.degrees(thetaDr)
    return thetaDd

# Calculate length of scattering vector
def Q_length(xlambda,thetaDd):
    Q = (2*PI/xlambda)*2*math.sin(math.radians(thetaDd))
    return Q

# Calculate theta using half-deflection-angle(thetaD), mu, gamma
def theta_angle(thetaDd,mud,gammad):
    thetaDr = math.radians(thetaDd)
    sinthetaD = math.sin(thetaDr)
    mur = math.radians(mud)
    gammar = math.radians(gammad)
    sinmu = math.sin(mur)
    cosmu = math.cos(mur)
    singamma = math.sin(gammar)
    cosgamma = math.cos(gammar)
    # When mu = 90 deg or gamma = 90 deg, thetaD is independent of theta
    if abs(cosmu*cosgamma) == 0:
        flag = False
        return flag, 0
    else:
        costwotheta = (1 - sinmu*singamma - 2*sinthetaD**2)/(cosmu*cosgamma)
    # Impossible for costwotheta > 1
    if abs(costwotheta) > 1:
        flag = False
        return flag, 0
    else:
        thetar = math.acos(costwotheta) / 2
        thetad = math.degrees(thetar)
        flag = True
        return flag, thetad

# Eta angle between hb (vector HKL) and nb (reference vector, generally surface normal)
# nb = np.array([[nH],[nK],[nL]], dtype=float)
def eta_angle(B,hb,nb):
    hc = np.dot(B,hb)
    nc = np.dot(B,nb)
    coseta = np.dot(np.transpose(hc), nc)[0][0] / (np.linalg.norm(hc)*np.linalg.norm(nc))
    etar = math.acos(coseta)
    etad = math.degrees(etar)
    return etad

# Get hb (vector HKL) using positions of angles
# H = hb[0][0], K = hb[1][0], L = hb[2][0]
# Note: theta = tth/2, not the position of motor th
def CheckBr(xlambda,U,B,thetad,omegad,chid,phid,mud,gammad):
    thetaDd = thetaD_angle_2(thetad,mud,gammad)
    sinthetaD = math.sin(math.radians(thetaDd))
    THETA_I = np.transpose(THETA_matrix(thetad))
    OMEGA_I = np.transpose(OMEGA_matrix(omegad))
    CHI_I = np.transpose(CHI_matrix(chid))
    PHI_I = np.transpose(PHI_matrix(phid))
    U_I = np.transpose(U)
    B_I = np.linalg.inv(B)
    D = D_matrix(thetad,mud,gammad)
    # Scattering vector in diffraction-plane space
    ksd = 2 * sinthetaD * (2*PI/xlambda) * np.array([[1], [0], [0]], dtype=float)
    hb = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(B_I,U_I),PHI_I),CHI_I),OMEGA_I),THETA_I),D),ksd)
    return hb

# Get azimuth angle psi, incident angle alphain, exit angle betaout
# nb: reference vector, generally surface normal
# nb = np.array([[nH],[nK],[nL]], dtype=float)
# Note: theta = tth/2, not the position of motor th
def CheckPsiAlphainBetaout(U,B,nb,thetad,omegad,chid,phid,mud,gammad):
   # Below to calculate psi
    THETA = THETA_matrix(thetad)
    OMEGA = OMEGA_matrix(omegad)
    CHI = CHI_matrix(chid)
    PHI = PHI_matrix(phid)
    D = D_matrix(thetad,mud,gammad)
    D_I = np.linalg.inv(D)
    # Reference vector in laboratory space
    nl = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(THETA,OMEGA),CHI),PHI),U),B),nb)
    # Reference vector in diffraction-plane space
    nd = np.dot(D_I,nl)
    # nd2, nd3 are y component and z component
    nd2 = nd[1][0]
    nd3 = nd[2][0]
    psir = math.atan2(nd3,nd2)
    psid = math.degrees(psir)
   # Below to calculate incident angle alphain
    mur = math.radians(mud)
    cosmu = math.cos(mur)
    sinmu = math.sin(mur)
    # mupl: unit vector of primary x-ray beam in laboratory space
    upl1 = 0
    upl2 = - cosmu
    upl3 = - sinmu
    mupl = (-1) * np.array([[upl1], [upl2], [upl3]], dtype=float)
    if np.linalg.norm(nl)*np.linalg.norm(mupl) != 0:
        sinalphain = np.dot(np.transpose(nl), mupl)[0][0] / (np.linalg.norm(nl)*np.linalg.norm(mupl))
    else:
        sinalphain = 0
    alphainr = math.asin(sinalphain)
    alphaind = math.degrees(alphainr)
    # Below to calculate exit angle betaout
    twothetar = math.radians(2*thetad)
    gammar = math.radians(gammad)
    sintwotheta = math.sin(twothetar)
    costwotheta = math.cos(twothetar)
    cosgamma = math.cos(gammar)
    singamma = math.sin(gammar)
    # udl: unit vector of deflected x-ray beam in laboratory space
    udl1 = sintwotheta * cosgamma
    udl2 = - costwotheta * cosgamma
    udl3 = - singamma
    udl = np.array([[udl1], [udl2], [udl3]], dtype=float)
    if np.linalg.norm(nl)*np.linalg.norm(udl) != 0:
        sinbetaout = np.dot(np.transpose(nl), udl)[0][0] / (np.linalg.norm(nl)*np.linalg.norm(udl))
    else:
        sinbetaout = 0
    betaoutr = math.asin(sinbetaout)
    betaoutd = math.degrees(betaoutr)
    return psid, alphaind, betaoutd

# Get the limitation of incident (exit) angles
# Note: theta = tth/2, not the position of motor th
# hb: vector HKL; nb: reference vector, generally surface normal
def CheckRangeAlphainBetaout(xlambda,B,hb,nb):
    mind = 0
    maxd = 0
    flag = False
    flagThetaD, thetaDd = thetaD_angle(xlambda,B,hb)
    if flagThetaD == False:
        return flag, mind, maxd
    etad = eta_angle(B,hb,nb)
    maxd = min(90, thetaDd + etad)
    mind = max(-90, thetaDd - etad)
    flag = True
    return flag, mind, maxd

# For a specified HKL, calculate four unknown angles from other three fixed angles
# Usage: e.g., to calculate th, theta, chi, phi from fixed mu, gamma, omega:
#     flag, angles = SC_angles(xlambda,U,B,hb,'x','x',omegad,'x','x',mud,gammad)
# flag: True (False) when calculation is successful (failed: no solution)
# angles: N set of positions, [[theta1, omega1, chi1, phi1, mu1, gamma1],[theta2,...],...,[thetaN,...]]
def SC_angles(xlambda,U,B,hb,thd,thetad,omegad,chid,phid,mud,gammad):
    angles = []
    unknown = []
    marker = []
    h = np.dot(np.dot(U,B),hb)
    h1 = h[0][0] / np.linalg.norm(h)
    h2 = h[1][0] / np.linalg.norm(h)
    h3 = h[2][0] / np.linalg.norm(h)
    flagThetaDd, thetaDd = thetaD_angle(xlambda,B,hb)
    if flagThetaDd == False:
        flag = False
        return flag, angles
    thetaDr = math.radians(thetaDd)
    costwothetaD = math.cos(thetaDr*2)
    sinthetaD = math.sin(thetaDr)
    udN = [thetad,mud,gammad].count('x')
    umN = [omegad,chid,phid].count('x')
    utN = [thd].count('x')
    uN = udN + umN
    # Four unknowns; at least one unknown in {theta, mu, gamma}
    if (udN>=1 and udN+umN+utN==4) == False:
        flag = False
        return flag, angles
    def equations(var):
        ui = -1
        if thetad == 'x':
            ui = ui+1
            thetar = math.radians(var[ui])
            marker.append(0)
        else:
            thetar = math.radians(thetad)
        sintheta = math.sin(thetar)
        costheta = math.cos(thetar)
        sintwotheta = math.sin(thetar*2)
        costwotheta = math.cos(thetar*2)
        if omegad == 'x':
            ui = ui+1
            omegar = math.radians(var[ui])
            marker.append(1)
        else:
            omegar = math.radians(omegad)
        sinomega = math.sin(omegar)
        cosomega = math.cos(omegar)
        if chid == 'x':
            ui = ui+1
            chir = math.radians(var[ui])
            marker.append(2)
        else:
            chir = math.radians(chid)
        sinchi = math.sin(chir)
        coschi = math.cos(chir)
        if phid == 'x':
            ui = ui+1
            phir = math.radians(var[ui])
            marker.append(3)
        else:
            phir = math.radians(phid)
        sinphi = math.sin(phir)
        cosphi = math.cos(phir)
        if mud == 'x':
            ui = ui+1
            mur = math.radians(var[ui])
            marker.append(4)
        else:
            mur = math.radians(mud)
        sinmu = math.sin(mur)
        cosmu = math.cos(mur)
        if gammad == 'x':
            ui = ui+1
            gammar = math.radians(var[ui])
            marker.append(5)
        else:
            gammar = math.radians(gammad)
        singamma = math.sin(gammar)
        cosgamma = math.cos(gammar)
        u1 = (sintwotheta*costheta*cosgamma + sintheta*cosmu - sintheta*costwotheta*cosgamma) / (2*sinthetaD)
        u2 = (-sintheta*sintwotheta*cosgamma + costheta*cosmu - costheta*costwotheta*cosgamma) / (2*sinthetaD)
        u3 = (sinmu - singamma) / (2*sinthetaD)
        R11 = cosphi*coschi*cosomega - sinphi*sinomega
        R12 = cosphi*coschi*sinomega + sinphi*cosomega
        R13 = cosphi*sinchi
        R21 = -sinphi*coschi*cosomega - cosphi*sinomega
        R22 = -sinphi*coschi*sinomega + cosphi*cosomega
        R23 = -sinphi*sinchi
        R31 = -sinchi*cosomega
        R32 = -sinchi*sinomega
        R33 = coschi
        eq1 = R11*u1 + R12*u2 + R13*u3 - h1
        eq2 = R21*u1 + R22*u2 + R23*u3 - h2
        eq3 = R31*u1 + R32*u2 + R33*u3 - h3
        eq4 = (1-sinmu*singamma-2*sinthetaD**2)/(cosmu*cosgamma) - costwotheta
        if thd == 'x':
            return [eq1,eq2,eq3,eq4]
        else:
            eq5 = math.radians(thd) - omegar - thetar
            return [eq1,eq2,eq3,eq4,eq5]
    solution = []
    for loop in range(0,180):
        ini = [(random.random()-0.5)*360 for i in range(0,uN)]
        sol = least_squares(equations,ini,bounds=(-180,180),gtol=1e-12)
        if sol.cost < 1e-6:
            solution.append(sol.x)
    # Remove duplicate results, precision 5 in output    
    solution = [[round(sij,5) for sij in si] for si in solution]
    solution = list(set([tuple(si) for si in solution]))
    if len(solution) == 0:
        flag = False
        return flag, angles
    angles = [[thetad,omegad,chid,phid,mud,gammad] for i in range(0,len(solution))]
    for i in range(0,len(solution)):
        for j in range(0,uN):
            angles[i][marker[j]] = solution[i][j]
    flag = True
    return flag, angles

# For a specified HKL and a desired psi (azimuth angle), calculate five unknown angles from other two fixed angles
# Usage: e.g., to calculate theta, omega, chi, phi from fixed mu, gamma:
#     flag, angles = SC_angles_fix_psi(xlambda,U,B,hb,nb,thd,psid,'x','x','x','x',mud,gammad)
# flag: True (False) when calculation is successful (failed: no solution)
# angles: N set of positions, [[theta1, omega1, chi1, phi1, mu1, gamma1],[theta2,...],...,[thetaN,...]]
def SC_angles_fix_psi(xlambda,U,B,hb,nb,psid,thd,thetad,omegad,chid,phid,mud,gammad):
    angles = []
    unknown = []
    marker = []
    hphi = np.dot(np.dot(U,B),hb)
    nphi = np.dot(np.dot(U,B),nb)
    t1phi = hphi
    t3phi = np.cross(t1phi.flatten(), nphi.flatten()).reshape(3,1)
    t2phi = np.cross(t3phi.flatten(), t1phi.flatten()).reshape(3,1)
    if np.linalg.norm(t3phi) == 0:
        flag = False
        return flag, angles
    t1phi = t1phi / np.linalg.norm(t1phi)
    t2phi = t2phi / np.linalg.norm(t2phi)
    t3phi = t3phi / np.linalg.norm(t3phi)
    Tphi = np.hstack((t1phi, t2phi, t3phi))
    Tphi_I = np.linalg.inv(Tphi)
    flagThetaDd, thetaDd = thetaD_angle(xlambda,B,hb)
    if flagThetaDd == False:
        flag = False
        return flag, angles
    thetaDr = math.radians(thetaDd)
    costwothetaD = math.cos(thetaDr*2)
    sinthetaD = math.sin(thetaDr)
    PSI = PSI_matrix(psid)
    udN = [thetad,mud,gammad].count('x')
    umN = [omegad,chid,phid].count('x')
    utN = [thd].count('x')
    uN = udN + umN
    if (udN+umN+utN==5) == False:
        flag = False
        return flag, angles
    def equations(var):
        ui = -1
        if thetad == 'x':
            ui = ui+1
            thetad_ = var[ui]
            marker.append(0)
        else:
            thetad_ = thetad
        if omegad == 'x':
            ui = ui+1
            omegad_ = var[ui]
            marker.append(1)
        else:
            omegad_ = omegad
        if chid == 'x':
            ui = ui+1
            chid_ = var[ui]
            marker.append(2)
        else:
            chid_ = chid
        if phid == 'x':
            ui = ui+1
            phid_ = var[ui]
            marker.append(3)
        else:
            phid_ = phid
        if mud == 'x':
            ui = ui+1
            mud_ = var[ui]
            marker.append(4)
        else:
            mud_ = mud
        if gammad == 'x':
            ui = ui+1
            gammad_ = var[ui]
            marker.append(5)
        else:
            gammad_ = gammad
        costwotheta = math.cos(math.radians(thetad_)*2)
        sinmu = math.sin(math.radians(mud_))
        cosmu = math.cos(math.radians(mud_))
        singamma = math.sin(math.radians(gammad_))
        cosgamma = math.cos(math.radians(gammad_))
        D = D_matrix(thetad_,mud_,gammad_)
        THETA = THETA_matrix(thetad_)
        THETA_I = np.transpose(THETA_matrix(thetad_))
        RIGHT = np.dot(np.dot(np.dot(THETA_I,D),PSI),Tphi_I)
        OMEGA = OMEGA_matrix(omegad_)
        CHI = CHI_matrix(chid_)
        PHI = PHI_matrix(phid_)
        R = np.dot(np.dot(OMEGA,CHI),PHI)
        eq1 = R[0][0] - RIGHT[0][0]
        eq2 = R[0][1] - RIGHT[0][1]
        eq3 = R[0][2] - RIGHT[0][2]
        eq4 = R[1][0] - RIGHT[1][0]
        eq5 = R[1][1] - RIGHT[1][1]
        eq6 = R[1][2] - RIGHT[1][2]
        eq7 = R[2][0] - RIGHT[2][0]
        eq8 = R[2][1] - RIGHT[2][1]
        eq9 = R[2][2] - RIGHT[2][2]
        eq10 = (1-sinmu*singamma-2*sinthetaD**2)/(cosmu*cosgamma) - costwotheta
        if thd == 'x':
            return [eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10]
        else:
            eq11 = math.radians(thd) - math.radians(omegad_) - math.radians(thetad_)
            return [eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11]
    solution = []
    for loop in range(0,180):
        ini = [(random.random()-0.5)*360 for i in range(0,uN)]
        sol = least_squares(equations,ini,bounds=(-180,180),gtol=1e-12)
        if sol.cost < 1e-6:
            solution.append(sol.x)
    # Remove duplicate results, precision 5 in output
    solution = [[round(sij,5) for sij in si] for si in solution]
    solution = list(set([tuple(si) for si in solution]))
    if len(solution) == 0:
        flag = False
        return flag, angles
    angles = [[thetad,omegad,chid,phid,mud,gammad] for i in range(0,len(solution))]
    for i in range(0,len(solution)):
        for j in range(0,uN):
            angles[i][marker[j]] = solution[i][j]
    flag = True
    return flag, angles

# For a specified HKL and a desired alphain (incident angle) or betaout (exit angle), calculate four unknown angles from other two fixed angles
# Usage: e.g., to calculate theta, omega, chi, phi from fixed betaout mu gamma:
#     flag, angles = SC_angles_fix_ab(xlambda,U,B,hb,nb,thd,'x',betaoutd,'x','x','x','x',mud,gammad)
# flag: True (False) when calculation is successful (failed: no solution)
# angles, N set of positions, [[theta1, omega1, chi1, phi1, mu1, gamma1],[theta2,...],...,[thetaN,...]]
def SC_angles_fix_ab(xlambda,U,B,hb,nb,alphaind,betaoutd,thd,thetad,omegad,chid,phid,mud,gammad):
    angles = []
    if [alphaind, betaoutd].count('x') != 1:
        flag = False
        return flag, angles
    flagThetaD, thetaDd = thetaD_angle(xlambda,B,hb)
    if flagThetaD == False:
        flag = False
        return flag, angles
    thetaDr = math.radians(thetaDd)
    sinthetaD = math.sin(thetaDr)
    costhetaD = math.cos(thetaDr)
    etad = eta_angle(B,hb,nb)
    etar = math.radians(etad)
    sineta = math.sin(etar)
    coseta = math.cos(etar)
    # Find corresponding two psi positions: psi1 and psi2
    if alphaind != 'x':
        alphainr = math.radians(alphaind)
        sinalphain = math.sin(alphainr)
        cospsi = (sinalphain - sinthetaD*coseta) / (costhetaD*sineta)
    if betaoutd != 'x':
        betaoutr = math.radians(betaoutd)
        sinbetaout = math.sin(betaoutr)
        cospsi = (sinthetaD*coseta - sinbetaout) / (costhetaD*sineta)
    if abs(cospsi) > 1:
        flag = False
        return flag, angles
    psi1r = math.atan2((1-cospsi**2)**0.5, cospsi)
    psi1d = math.degrees(psi1r)
    psi2r = math.atan2(-(1-cospsi**2)**0.5, cospsi)
    psi2d = math.degrees(psi2r)
    flag1, angles1 = SC_angles_fix_psi(xlambda,U,B,hb,nb,psi1d,thd,thetad,omegad,chid,phid,mud,gammad)
    flag2, angles2 = SC_angles_fix_psi(xlambda,U,B,hb,nb,psi2d,thd,thetad,omegad,chid,phid,mud,gammad)
    angles = angles1 + angles2
    # Remove duplicate positions for psi1 and psi2
    angles = list(set([tuple(ai) for ai in angles]))
    if len(angles) == 0:
        flag = False
    else:
        flag = True
    return flag, angles

# Change (theta,omega,chi,phid,mu,gamma) to (tth,th,chi,phi,mu,gamma)
# Display angle positions in [-180,180)
# Sort N sets of positions
# tth = theta*2
# th = tth/2 + omega
def motors(angles):
    for i in range(0,len(angles)):
        thetad,omegad,chid,phid,mud,gammad = angles[i]
        tthd = thetad * 2
        thd = thetad + omegad
        angles[i] = [tthd,thd,chid,phid,mud,gammad]
        for j in range(0,6):
            while angles[i][j] < -180:
                angles[i][j] = angles[i][j] + 360.0
            while angles[i][j] >= 180:
                angles[i][j] = angles[i][j] - 360.0
            angles[i][j] = float(round(angles[i][j],5))
    # Remove duplicate results
    angles = list(set([tuple(ai) for ai in angles]))
    N = len(angles)
    # Sort N sets of positions, in the order of tth, th, chi, phi, mu, gamma
    index = [0 for i in range(0,N)]
    for i in range(0,N):
        for j in range(0,6):
            index[i] = index[i] + int(abs(angles[i][j])/18)*10**(5-j)
    angles =[a for i,a in sorted(zip(index,angles))]
    return N, angles
