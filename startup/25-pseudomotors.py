from ophyd import (EpicsMotor, PseudoSingle, PseudoPositioner, Device, Signal,
                   Component as Cpt)
from ophyd.pseudopos import (pseudo_position_argument, real_position_argument)
import numpy as np
from scipy import interpolate
from ophyd import SoftPositioner

from utils.sixcircle import SixCircle
import time

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
    

def hkl_to_angles(H, K, L):
    flag, pos = sc.ca_s(H, K, L)
    if flag == True:
        tth, th, chi, phi, caMU, caGAM, caSA, caOMEGA, caAZIMUTH, caALPHA, caBETA = pos[0]
        # sc.br(H, K, L)
        return tth, th, chi, phi
    else:
        return None
    

def angles_to_hkl(Tth, Th, Chi, Phi):
    H, K, L = sc.mv(tth=Tth,th=Th,chi=Chi,phi=Phi)
    return H, K, L


# class HKLPseudo(PseudoPositioner):
# # Defines H, K, L pseudomotors for the Six Circle diffractometer code
#     H = Cpt(PseudoSingle)
#     K = Cpt(PseudoSingle)
#     L = Cpt(PseudoSingle)

#     th = Cpt(SoftPositioner, kind="hinted", init_pos=0)
#     chi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
#     phi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
#     tth = Cpt(SoftPositioner, kind="hinted", init_pos=0)

    # th = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:Th}Mtr', labels=('scir',))
    # chi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:ChiA}Mtr', labels=('scir',))
    # phi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:PhiA}Mtr', labels=('scir',))
    # tth = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:2Th}Mtr', labels=('scir',))

    # def __init__(self, *, calc_to_real, real_to_calc, **kwargs):
    #     """
    #     calc_to_real: function (H, K, L) -> (tth, th, chi, phi)
    #     real_to_calc: function (tth, th, chi, phi) -> (H, K, L)
    #     """
    #     self._calc_to_real = calc_to_real
    #     self._real_to_calc = real_to_calc
    #     super().__init__(**kwargs)
    # @property
    # def position(self):
    #     """Return the current pseudo position (H, K, L)."""
    #     realpos = self.real_position
    #     h, k, l = angles_to_hkl(realpos.tth, realpos.th, realpos.chi, realpos.phi)
    #     sc.wh_refresh()
    #     return self.PseudoPosition(H=h, K=k, L=l)

    # @pseudo_position_argument
    # def forward(self, pseudopos):
    #     """(H, K, L) -> (tth, th, chi, phi)"""
    #     # print("H, K, L in pseduopos", H, K, L)
    #     res = hkl_to_angles(pseudopos.H, pseudopos.K, pseudopos.L)

    #     if res is None:
    #         raise ValueError(f"Cannot reach (H,K,L)=({pseudopos.h},{pseudopos.k},{pseudopos.l})")
        
    #     catth, cath, cachi, caphi = res
    #     return self.RealPosition(tth=catth, th=cath, chi=cachi, phi=caphi)
    
    # @real_position_argument
    # def inverse(self, realpos):
    #     """(tth, th, chi, phi) -> (H, K, L)"""
    #     # tth, th, chi, phi = realpos.Tth, realpos.Th, realpos.Chi, realpos.Phi
    #     # print("angles in realpos", tth, th, chi, phi)
    #     caH, caK, caL = angles_to_hkl(realpos.Tth, realpos.Th, realpos.Chi, realpos.Phi)
    #     return self.PseudoPosition(H=caH, K=caK, L=caL)


class HKLPseudoV2(PseudoPositioner):
    H = Cpt(PseudoSingle)
    K = Cpt(PseudoSingle)
    L = Cpt(PseudoSingle)

    th = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:Th}Mtr', labels=('scir',))
    chi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:ChiA}Mtr', labels=('scir',))
    phi = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:PhiA}Mtr', labels=('scir',))
    tth = Cpt(EpicsMotor, 'XF:10IDD-OP{Spec:1-Ax:2Th}Mtr', labels=('scir',))

    # th = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    # chi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    # phi = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    # tth = Cpt(SoftPositioner, kind="hinted", init_pos=0)
    
    def __init__(self, *args, **kwargs):
        self.sc = SixCircle()       
        super().__init__(**kwargs)

    def hkl_to_angles(self, H, K, L):
        # Optional: reject (0,0,0) early if invalid
        if np.allclose([H, K, L], [0, 0, 0]):
            return None  # or raise?
        
        flag, pos = self.sc.ca_s(H, K, L)
        if flag == True:
            tth, th, chi, phi, caMU, caGAM, caSA, caOMEGA, caAZIMUTH, caALPHA, caBETA = pos[0]
            # sc.br(H, K, L)
            return tth, th, chi, phi
        else:
            return None

    def angles_to_hkl(self, Tth, Th, Chi, Phi):
        # breakpoint()
        H, K, L = self.sc.mv(tth=Tth, th=Th, chi=Chi, phi=Phi)
        return H, K, L

    @property
    def position(self):
        """Return the current pseudo position (H, K, L)."""
        realpos = self.real_position
        h, k, l = self.angles_to_hkl(realpos.tth, realpos.th, realpos.chi, realpos.phi)
        self.sc.wh_refresh()
        return self.PseudoPosition(H=h, K=k, L=l)

    @pseudo_position_argument
    def forward(self, pseudopos):
        """(H, K, L) -> (tth, th, chi, phi)"""
        # print("H, K, L in pseduopos", H, K, L)
        res = self.hkl_to_angles(pseudopos.H, pseudopos.K, pseudopos.L)
        self.sc.wh_refresh()
        if res is None:
            raise ValueError(f"Cannot reach (H,K,L)=({pseudopos.h},{pseudopos.k},{pseudopos.l})")
        
        catth, cath, cachi, caphi = res
        return self.RealPosition(tth=catth, th=cath, chi=cachi, phi=caphi)
    
    @real_position_argument
    def inverse(self, realpos):
        """(tth, th, chi, phi) -> (H, K, L)"""
        # tth, th, chi, phi = realpos.Tth, realpos.Th, realpos.Chi, realpos.Phi
        # print("angles in realpos", tth, th, chi, phi)
        caH, caK, caL = self.angles_to_hkl(realpos.tth, realpos.th, realpos.chi, realpos.phi)
        self.sc.wh_refresh()
        return self.PseudoPosition(H=caH, K=caK, L=caL)
    

class HKLDerived(Device):
    alpha = Cpt(Signal, value=0.0, kind='hinted')
    beta  = Cpt(Signal, value=0.0, kind='hinted')
    omega = Cpt(Signal, value=0.0, kind='hinted')
    mu = Cpt(Signal, value=0.0, kind='hinted')
    gam = Cpt(Signal, value=0.0, kind='hinted')
    sa = Cpt(Signal, value=0.0, kind='hinted')
    azimuth = Cpt(Signal, value=0.0, kind='hinted')

    def __init__(self, hkl_pseudo, *args, **kwargs):
        # self.hkl_pseudo = hkl_pseudo
        super().__init__(*args, **kwargs)

        # Subscribe to changes in key real motors
        # for m in [self.hkl_pseudo.tth, self.hkl_pseudo.th,
        #           self.hkl_pseudo.chi, self.hkl_pseudo.phi]:
        #     m.subscribe(self._on_motor_change)

        # Initialize derived values
        self._update_from_motors(hkl_pseudo)

    def _on_motor_change(self, *args, **kwargs):
        """Callback when any subscribed motor position changes."""
        self._update_from_motors()

    def _update_from_motors(self, hkl_pseudo):
        """Recalculate derived parameters based on current positions."""
        try:
            hkl_pseudo = getattr(self, "hkl_pseudo", None)
            if hkl_pseudo is None:
                print("[HKLDerived] No hkl_pseudo device attached — skipping update.")
                return

            realpos = hkl_pseudo.real_position

            # Extract positions safely
            tth = getattr(realpos, "tth", None)
            th  = getattr(realpos, "th", None)
            chi = getattr(realpos, "chi", None)
            phi = getattr(realpos, "phi", None)

            if any(v is None for v in [tth, th, chi, phi]):
                print("[HKLDerived] Skipping update: one or more motor positions are None.")
                return

            if any(np.isnan(v) for v in [tth, th, chi, phi]):
                print("[HKLDerived] Skipping update: one or more motor positions are NaN.")
                return

            H, K, L = hkl_pseudo.angles_to_hkl(tth, th, chi, phi)

            if (H is None or K is None or L is None or
                np.isnan(H) or np.isnan(K) or np.isnan(L) or
                np.allclose([H, K, L], [0, 0, 0])):
                print("[HKLDerived] Skipping update: invalid or undefined HKL position.")
                return

            flag, pos = hkl_pseudo.sc.ca_s(H, K, L)
            if not flag or not pos:
                print("[HKLDerived] Skipping update: ca_s() returned invalid result.")
                return

            *_, caMU, caGAM, caSA, caOMEGA, caAZIMUTH, caALPHA, caBETA = pos[0]

            self.alpha.put(caALPHA)
            self.beta.put(caBETA)
            self.omega.put(caOMEGA)
            self.mu.put(caMU)
            self.gam.put(caGAM)
            self.sa.put(caSA)
            self.azimuth.put(caAZIMUTH)

        except Exception as ex:
            print(f"[HKLDerived] Update failed: {ex}")

    def trigger(self):
        """Refresh before any scan baseline read."""
        self._update_from_motors()
        return super().trigger()


anc_xtal = AnalyzerCXtal('', name='anc_xtal', egu=('deg', 'mm'))
sam_prime = SamplePrime('', name='sp', egu=('mm', 'deg'))
# hklpseudo = HKLPseudo(name='hkl', calc_to_real=hkl_to_angles, real_to_calc=angles_to_hkl)
hklps = HKLPseudoV2(name='hkl')
for m in [hklps.tth, hklps.th, hklps.chi, hklps.phi]:
    try:
        m.wait_for_connection(timeout=5)
        m.read()  # Force PV read to populate .position
    except Exception as ex:
        print(f"[Startup] Warning: {m.name} not ready ({ex})")

hkl_params = HKLDerived(hklps, name='hkl_params')
# hkl_params._update_from_motors(hklps)  # initial update