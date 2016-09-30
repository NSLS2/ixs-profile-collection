from ophyd import (EpicsMotor, PseudoSingle, PseudoPositioner,
                   Component as Cpt)
from ophyd.pseudopos import (pseudo_position_argument, real_position_argument)
import numpy as np
from scipy import interpolate


_ivu_gap = [ 6184.0,
       6368.0,
       6550.0,
       6730.0,
       6909.0,
       7086.0,
       7262.0,
       7438.0,
       7613.0,
       7788.0,
       7962.0,
       8138.0,
       8311.0,
       8486.0,
       8662.0,
       8840.0,
       9019.0,
       9200.0,
       9382.0,
       9567.0,
       9753.0,
       9943.0,
       10136.0,
       10332.0,
       10531.0,
       10740.0,
       10960.0,
       11385.0]

_Bragg = [ 16.401781,
       15.821071,
       15.280587,
       14.776273,
       14.304428,
       13.862182,
       13.446656,
       13.055574,
       12.686789,
       12.338436,
       12.008851,
       11.696607,
       11.400018,
       11.118472,
       10.850483,
       10.595196,
       10.351707,
       10.119222,
       9.896987,
       9.684357,
       9.480717,
       9.28554,
       9.098203,
       8.918375,
       8.7455,
       8.579178,
       8.41915,
       8.11645 ]

gcalc = interpolate.interp1d(_Bragg, _ivu_gap)

_hc = 12398.4193
_si_111 = 3.1363


class BLEnergy(PseudoPositioner):
    # limits and constants from Spec file, "site.mac"
    energy = Cpt(PseudoSingle, limits=(7.835, 17.7))

    theta = Cpt(EpicsMotor, 'XF:10IDA-OP{Mono:DCM-Ax:P}Mtr')
    z2 = Cpt(EpicsMotor, 'XF:10IDA-OP{Mono:DCM-Ax:Z2}Mtr')
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

    theta = Cpt(EpicsMotor, 'XF:10IDA-OP{Mono:DCM-Ax:P}Mtr')
    z2 = Cpt(EpicsMotor, 'XF:10IDA-OP{Mono:DCM-Ax:Z2}Mtr')


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

    uof = Cpt(EpicsMotor, 'XF:10IDB-OP{Mono:HRM2-Ax:UTO}Mtr')
    dof = Cpt(EpicsMotor, 'XF:10IDB-OP{Mono:HRM2-Ax:DTO}Mtr')

    def forward(self, pseudo_pos):
        _pos = -1.0*pseudo_pos.energy*np.tan(np.deg2rad(_TB))/_EB
        return self.RealPosition(uof=_pos, dof=_pos)

    def inverse(self, real_pos):
        _energy = -1.0*_EB*real_pos.uof/np.tan(np.deg2rad(_TB))
        return self.PseudoPosition(energy=_energy)


hrmE = HRMEnergy('', name='hrmE', egu='eV')


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
b2 = bny*bny + bnz*bnz
cny = 2.72030
cnz = 3.40940
c2 = cny*cny + cnz*cnz
c1 = np.sqrt(c2)
d2 = B0y*B0y + B0z*B0z
dn = np.sqrt(d2)
h2 = C0y*C0y + C0z*C0z
hn = np.sqrt(h2)


class AnalyzerCXtal(PseudoPositioner):
    the  = Cpt(PseudoSingle, egu='deg')
    y = Cpt(PseudoSingle, egu='mm')

    uy = Cpt(EpicsMotor, 'XF:10IDD-OP{Analy:1-Ax:UY}Mtr')
    dy = Cpt(EpicsMotor, 'XF:10IDD-OP{Analy:1-Ax:DY}Mtr')

    @pseudo_position_argument
    def forward(self, pseudopos):
        ddy = np.zeros((2,), dtype=float)
        cth = pseudopos.the + th0
        ccy = pseudopos.y

        cy = c1*np.sin(np.deg2rad(cth))
        cz = np.sqrt(c2 - cy*cy)
        a1 = (c2 + h2 - d2)/2.
        d1 = cz*np.sqrt(c2*h2 - a1*a1)
        hy = (a1*cy + d1)/c2
        #
        # first motor position (mm)
        y1 = 0.01*ccy - hy - C0y
        ddy[0] = 100.*y1
        #
        a2 = anz*anz + c2 - b2 - 2*anz*cz
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
        d1y = self.uy.position
        d2y = self.dy.position

        _any = 0.84992 + 0.01*(d2y - d1y)
        a2 = _any*_any + anz*anz
        a1 = (c2 + a2 - b2)/2.
        d1 = anz*np.sqrt(a2*c2 - a1*a1)
        cy = (a1*_any + d1)/a2
        cz = np.sqrt(c2 - cy*cy)
        #
        # crystal angle (rad)
        ccp[0] = np.rad2deg(np.arcsin(cy/c1))
        a2 = (c2 + h2 - d2)/2.
        _d2 = cz*np.sqrt(c2*h2 - a2*a2)
        hy = (a2*cy + _d2)/c2
        #
        # crystal y-position
        ccp[1] = 100.*(C0y + 0.01*d1y + hy)

        return self.PseudoPosition(the=ccp[0]-th0, y=ccp[1])


anc_xtal = AnalyzerCXtal('', name='anc_xtal', egu=('deg', 'mm'))
