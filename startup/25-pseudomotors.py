from ophyd import (EpicsMotor, PseudoSingle, PseudoPositioner,
                   Component as Cpt)
from ophyd.pseudopos import (pseudo_position_argument, real_position_argument)
import numpy as np
from scipy import interpolate

# List of available EpicsMotor labels in this script
# [blenergy, dcmenergy, hrmenergy, analyzercxtal]

_ivu_gap = [6184.0, 6368.0, 6550.0, 6730.0, 6909.0, 7086.0, 7262.0,
            7438.0, 7613.0, 7788.0, 7962.0, 8138.0, 8311.0, 8486.0,
            8662.0, 8840.0, 9019.0, 9200.0, 9382.0, 9567.0, 9753.0,
            9943.0, 10136.0, 10332.0, 10531.0, 10740.0, 10960.0,
            11385.0]

_Bragg = [16.401781, 15.821071, 15.280587, 14.776273, 14.304428,
          13.862182, 13.446656, 13.055574, 12.686789, 12.338436,
          12.008851, 11.696607, 11.400018, 11.118472, 10.850483,
          10.595196, 10.351707, 10.119222, 9.896987, 9.684357, 9.480717,
          9.28554, 9.098203, 8.918375, 8.7455, 8.579178, 8.41915, 8.11645]

gcalc = interpolate.interp1d(_Bragg, _ivu_gap)

_hc = 12398.4193
_si_111 = 3.1363


class BLEnergy(PseudoPositioner):
    # limits and constants from Spec file, "site.mac"
    energy = Cpt(PseudoSingle, limits=(7.835, 17.7))

    theta = Cpt(EpicsMotor, 'XF:10IDA-OP{Mono:DCM-Ax:P}Mtr', labels=('blenergy',))
    z2 = Cpt(EpicsMotor, 'XF:10IDA-OP{Mono:DCM-Ax:Z2}Mtr', labels=('blenergy',))
    ugap = Cpt(Undulator, 'SR:C10-ID:G1{IVU22:1')

    def forward(self, pseudo_pos):
        _th = np.rad2deg(np.arcsin(_hc/(2.*_si_111*pseudo_pos)))
        _z2 = 15./np.cos(np.deg2rad(_th)) - 15.
        _ugap = gcalc(_th)

        return self.RealPosition(theta=_th, z2=_z2, ugap=_ugap)

    def inverse(self, real_pos):
        en = _hc/(2.*_si_111*np.sin(np.deg2rad(real_pos.theta)))
        en = float(en)
        return self.PseudoPosition(energy=en)


blE = BLEnergy('', name='blE', egu='eV')


class DCMEnergy(PseudoPositioner):
    # limits and constants from Spec file, "site.mac"
    energy = Cpt(PseudoSingle, limits=(7.835, 17.7))

    theta = Cpt(EpicsMotor, 'XF:10IDA-OP{Mono:DCM-Ax:P}Mtr', labels=('dcmenergy',))
    z2 = Cpt(EpicsMotor, 'XF:10IDA-OP{Mono:DCM-Ax:Z2}Mtr', labels=('dcmenergy',))

    def forward(self, pseudo_pos):
        _th = np.rad2deg(np.arcsin(_hc/(2.*_si_111*pseudo_pos)))
        _z2 = 15./np.cos(np.deg2rad(_th)) - 15.

        return self.RealPosition(theta=_th, z2=_z2)

    def inverse(self, real_pos):
        en = _hc/(2.*_si_111*np.sin(np.deg2rad(real_pos.theta)))
        en = float(en)
        return self.PseudoPosition(energy=en)


dcmE = DCMEnergy('', name='dcmE', egu='eV')


# constants for HRM2 energy conversion
_TB = 69.2943
_EB = 9.1317


class HRMEnergy(PseudoPositioner):
    energy = Cpt(PseudoSingle)

    uof = Cpt(EpicsMotor, 'XF:10IDB-OP{Mono:HRM2-Ax:UTO}Mtr', labels=('hrmenergy',))
    dof = Cpt(EpicsMotor, 'XF:10IDB-OP{Mono:HRM2-Ax:DTO}Mtr', labels=('hrmenergy',))

    def forward(self, pseudo_pos):
        _pos = -1.0*pseudo_pos.energy*np.tan(np.deg2rad(_TB))/_EB
        return self.RealPosition(uof=_pos, dof=_pos)

    def inverse(self, real_pos):
        _energy = -1.0*_EB*real_pos.uof/np.tan(np.deg2rad(_TB))
        return self.PseudoPosition(energy=_energy)


hrmE = HRMEnergy('', name='hrmE', egu='eV')
hrmE.energy.readback.name = hrmE.name
hrmE.uof.user_readback.kind = 'normal'
hrmE.dof.user_readback.kind = 'normal'

# ripped from Spec file, "site.mac"
th0 = 38.5857476

A0y = 1.21523
A0z = 7.20726
B0y = 3.08561
B0z = 6.57882
C0y = 0.36531
C0z = 3.16942
anz = 4.03784
bny = 1.87038
bnz = -0.62844
b2_par = bny*bny + bnz*bnz
cny = 2.72030
cnz = 3.40940
c2_par = cny*cny + cnz*cnz
c1_par = np.sqrt(c2_par)
d2_par = B0y*B0y + B0z*B0z
dn = np.sqrt(d2_par)
h2_par = C0y*C0y + C0z*C0z
hn_par = np.sqrt(h2_par)


class AnalyzerCXtal(PseudoPositioner):
    the = Cpt(PseudoSingle, egu='deg')
    y = Cpt(PseudoSingle, egu='mm')

    uy = Cpt(EpicsMotor, 'XF:10IDD-OP{Analy:1-Ax:UY}Mtr', labels=('analyzercxtal',))
    dy = Cpt(EpicsMotor, 'XF:10IDD-OP{Analy:1-Ax:DY}Mtr', labels=('analyzercxtal',))

    @pseudo_position_argument
    def forward(self, pseudopos):
        ddy = np.zeros((2,), dtype=float)
        cth = pseudopos.the + th0
        ccy = pseudopos.y

        cy = c1_par*np.sin(np.deg2rad(cth))
        cz = np.sqrt(c2_par - cy*cy)
        a1 = (c2_par + h2_par - d2_par)/2.
        d1 = cz*np.sqrt(c2_par*h2_par - a1*a1)
        hy = (a1*cy + d1)/c2_par
        #
        # first motor position (mm)
        y1 = 0.01*ccy - hy - C0y
        ddy[0] = 100.*y1
        #
        a2 = anz*anz + c2_par - b2_par - 2*anz*cz
        d3 = np.sqrt(cy*cy - a2)
        #
        # second motor position (mm)
        y2 = C0y + y1 + cy - d3 - A0y
        ddy[1] = 100.*y2

        return self.RealPosition(uy=ddy[0], dy=ddy[1])

    @real_position_argument
    def inverse(self, realpos):
        ccp = np.zeros((2,), dtype=float)
        # d2y == a2uy, d1y == a2dy
        #TODO use realpos input!
        d1y = self.uy.position
        d2y = self.dy.position

        _any = 0.84992 + 0.01*(d2y - d1y)
        a2 = _any*_any + anz*anz
        a1 = (c2_par + a2 - b2_par)/2.
        d1 = anz*np.sqrt(a2*c2_par - a1*a1)
        cy = (a1*_any + d1)/a2
        cz = np.sqrt(c2_par - cy*cy)
        #
        # crystal angle (rad)
        ccp[0] = np.rad2deg(np.arcsin(cy/c1_par))
        a2 = (c2_par + h2_par - d2_par)/2.
        _d2 = cz*np.sqrt(c2_par*h2_par - a2*a2)
        hy = (a2*cy + _d2)/c2_par
        #
        # crystal y-position
        ccp[1] = 100.*(C0y + 0.01*d1y + hy)

        return self.PseudoPosition(the=ccp[0]-th0, y=ccp[1])
    

class SamplePrime(PseudoPositioner):
    xp = Cpt(PseudoSingle, egu='mm')
    zp = Cpt(PseudoSingle, egu='mm')

    th = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:Th}Mtr', labels=('sampleprime',))
    phi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:PhiA}Mtr', labels=('sampleprime',))
    sx = Cpt(EpicsMotor, 'XF:10IDD-OP{Env:1-Ax:X}Mtr', labels=('sampleprime',))
    sz = Cpt(EpicsMotor, 'XF:10IDD-OP{Env:1-Ax:Z}Mtr', labels=('sampleprime',))

    @pseudo_position_argument
    def forward(self, pseudopos):
        _th = np.deg2rad(self.th.position + self.phi.position)
        _sx = pseudopos.xp*np.cos(_th) - pseudopos.zp*np.sin(_th)
        _sz = pseudopos.xp*np.sin(_th) + pseudopos.zp*np.cos(_th)

        return self.RealPosition(sx = _sx, sz = _sz)
    
    @real_position_argument
    def inverse(self, realpos):
        _th = np.deg2rad(realpos.th + realpos.phi)
        _xp = realpos.sz*np.sin(_th) + realpos.sx*np.cos(_th)
        _zp = realpos.sz*np.cos(_th) - realpos.sx*np.sin(_th)

        return self.PseudoPosition(xp = _xp, zp = _zp)
    

# class HKLPseudo(PseudoPositioner):
# # Defines H, K, L pseudomotors for the Six Circle diffractometer code
#     H = Cpt(PseudoSingle)
#     K = Cpt(PseudoSingle)
#     L = Cpt(PseudoSingle)

#     th = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:Th}Mtr', labels=('scir',))
#     chi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:ChiA}Mtr', labels=('scir',))
#     phi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:PhiA}Mtr', labels=('scir',))
#     tth = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:2Th}Mtr', labels=('scir',))

#     def __init__(self, *, calc_to_real, real_to_calc, **kwargs):
#         """
#         calc_to_real: function (H, K, L) -> (tth, th, chi, phi)
#         real_to_calc: function (tth, th, chi, phi) -> (H, K, L)
#         """
#         self._calc_to_real = calc_to_real
#         self._real_to_calc = real_to_calc
#         super().__init__(**kwargs)

#     @pseudo_position_argument
#     def forward(self, pseudopos):
#         """(H, K, L) -> (tth, th, chi, phi)"""
#         H, K, L = pseudopos
#         catth, cath, cachi, caphi = self._calc_to_real(H, K, L)
#         return self.RealPosition(tth=catth, th=cath, chi=cachi, phi=caphi)
    
#     @real_position_argument
#     def inverse(self, realpos):
#         """(tth, th, chi, phi) -> (H, K, L)"""
#         tth, th, chi, phi = realpos
#         caH, caK, caL = self._real_to_calc(tth, th, chi, phi)
#         return self.PseudoPosition(H=caH, K=caK, L=caL)


def hkl_to_angles(H, K, L):
    flag, pos = sc.ca_s(H, K, L)
    if flag == True:
        tth, th, chi, phi, caMU, caGAM, caSA, caOMEGA, caAZIMUTH, caALPHA, caBETA = pos[0]
        return tth, th, chi, phi
    else:
        return
    

def angles_to_hkl(Tth, Th, Chi, Phi):
    sc.mv(tth=Tth,th=Th,chi=Chi,phi=Phi)
    H, K, L = sc.wh_refresh()
    return H, K, L



anc_xtal = AnalyzerCXtal('', name='anc_xtal', egu=('deg', 'mm'))
sam_prime = SamplePrime('', name='sp', egu=('mm', 'deg'))
# hkl = HKLPseudo(name='hkl', calc_to_real=hkl_to_angles, real_to_calc=angles_to_hkl)