import bluesky.plans as bp
import bluesky.preprocessors as bpp
import bluesky.plan_stubs as bps
from bluesky.callbacks.fitting import PeakStats


#*******************************************************************************************************
def qq2th(Qq,En=9.1317):
    # Calculates scattering TH-angle from the Q-value. Energy En is set in keV.
    if En > 15.0:
       print("Wrong energy value")
       return

    Wl = 1.24/En
    Th = np.degrees(2.0*np.arcsin(Wl*Qq/(4.0*np.pi)))
    return Th


#*******************************************************************************************************
def th2qq(Th,En=9.1317):
    # Calculates the Q-value from the scattering TH-angle. Energy En is set in keV.
    if En > 15.0:
       print("Wrong energy value")
       return

    Wl = 1.24/En
    Qq = 4.*np.pi*np.sin(np.radians(Th/2.))/Wl
    return Qq


#def set_acquire_time(t):
#    yield from bps.mv(sclr.preset_time, t)
#    yield from bps.mv(lambda_det.cam.acquire_time, t*0.995)
#    yield from bps.mv(lambda_det.cam.acquire_period, t*0.995)


#*******************************************************************************************************
def hrmE_dscan(start, stop, steps, exp_time, md=None):
    """
    Run a relative (delta) hmre scan with lambda and scalar

    Paramameters
    ------------
    start, stop : float
        The relative start and stop points

    steps : int
        The number of gaps (take steps + 1 measurements)

    t : float
        The exposure time in seconds
    """

    md = md or {}
    md['count_time'] = exp_time

    yield from set_lambda_exposure(exp_time)
    return (
#        yield from dscan(hrmE, start, stop, steps, [lambda_det], det_channel=[6])
        yield from bp.rel_scan([lambda_det], hrmE, start, stop, steps, md=md)
    )


#*******************************************************************************************************
def hrmE_ascan(start, stop, steps, exp_time, md=None):
    """
    Run a absolute hmre scan with lambda and scalar

    Paramameters
    ------------
    start, stop : float
        The absolute start and stop points

    steps : int
        The number of gaps (take steps + 1 measurements)

    t : float
        The exposure time in seconds
    """

    md = md or {}
    md['count_time'] = exp_time

    yield from set_lambda_exposure(exp_time)
    return (
        yield from bp.scan([lambda_det], hrmE, start, stop, steps, md=md)
    )

#*******************************************************************************************************
from ophyd.sim import SynAxis, SynGauss

# Simulated motor (positioner)
sym_mot = SynAxis(name="sym_mot")

# Simulated detector with Gaussian response vs motor position
#    I(m) = Imax * exp(-(m - center)^2 / (2*sigma^2)) + noise
sym_det = SynGauss(
    name="sym_det",
    motor=sym_mot,
    motor_field="sym_mot",   # must match key in sym_mot.read()
    center=0.0,
    Imax=1e5,
    sigma=0.2,
    noise="poisson",
    random_state=0,
)
