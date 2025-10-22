'''
from bluesky.callbacks.olog import OlogCallback
from bluesky.callbacks import LivePlot


from suitcase.spec import DocumentToSpec

spec_cb = DocumentToSpec('/tmp/spec1.spec')
RE.subscribe(spec_cb)
'''
import tkinter as tk
from tkinter import filedialog

from utils.sixcircle import SixCircle

sc = SixCircle()

# def whsc():
# # Defines the current position of the sample
#     Th = spec.th.position
#     Tth = spec.tth.position
#     Chi = spec.chi.position
#     Phi = spec.phi.position
#     sc.mv(tth=Tth,th=Th,chi=Chi,phi=Phi)


def wh():
    # updates the H, K, L pseudomotors from the current position of the sample
    hklps.sc.wh()


# def wh_refresh():
#     res = hklps.position
#     sc.wh_refresh()
#     # return sc.H, sc.K, sc.L


def freeze(*arg):
    hklps.sc.freeze(*arg)


def setfrozen(*arg):
    hklps.sc.setfrozen(*arg)


def pa():
    hklps.sc.pa()


def or0(h,k,l):
# defines the OR0 vector from the current angle positions
    hklps.sc.or0(h,k,l)


def or1(h,k,l):
# defines the OR1 vector from the current angle positions
    hklps.sc.or1(h,k,l)


# def update_spec_signals():
#    spec.H.put(sc.H)
#    spec.K.put(sc.K)
#    spec.L.put(sc.L)
#    spec.Q.put(sc.ABSQ*10)
#    spec.LAMBDA.put(sc.LAMBDA)
#    spec.HAZ.put(sc.g_haz)
#    spec.KAZ.put(sc.g_kaz)
#    spec.LAZ.put(sc.g_laz)
#    spec.AZIMUTH.put(sc.AZIMUTH)
#    spec.ALPHA.put(sc.ALPHA)
#    spec.BETA.put(sc.BETA)
#    spec.OMEGA.put(sc.OMEGA)


# def br(*args):
#     hklps.sc.br(*args)
    # update_spec_signals()
    #yield from bps.mv(spec.th, sc.TH, spec.tth, sc.TTH, spec.chi, sc.CHI, spec.phi, sc.PHI)

def br(h, k, l):
    # moves the sample to the desired H, K, L positions
    yield from bps.mv(hklps.H, h, hklps.K, k, hklps.L, l)


def mv(mu=None, gam=None, tth=None, th=None, chi=None, phi=None):
    """
    Smart move dispatcher:
      - For mu, gam: move via SixCircle (geometry engine).
      - For tth, th, chi, phi: move real motors using Bluesky RunEngine.
    """
    moves = {}  # collect real motor moves
    sc_args = {}  # for self.sc.mv()

    # classify arguments
    if mu is not None:
        sc_args['mu'] = mu
    if gam is not None:
        sc_args['gam'] = gam
    if tth is not None:
        moves[hklps.tth] = tth
    if th is not None:
        moves[hklps.th] = th
    if chi is not None:
        moves[hklps.chi] = chi
    if phi is not None:
        moves[hklps.phi] = phi

    # move motors first (if any)
    if moves:
        yield from bps.mv(*sum(moves.items(), ()))  # expands {motor: pos} → motor, pos, motor, pos
    # move the geometry (mu/gam) next
    if sc_args:
        hklps.sc.mv(**sc_args)


def ca(h, k, l):
    res = hklps.sc.ca(h,k,l)


def select_file_to_open(initialdir=".", filetypes=(("All files", "*.*"),)):
    # Opens a file dialog to select a file
    root = tk.Tk()
    root.withdraw()  # hide the empty Tk window
    filename = filedialog.askopenfilename(initialdir=initialdir, filetypes=filetypes)
    root.destroy()
    return filename

def select_file_to_save(initialdir="."):
    # Opens a file dialog to select a directory
    root = tk.Tk()
    root.withdraw()
    filename = filedialog.asksaveasfilename(initialdir=initialdir, defaultextension=".conf", filetypes=(("Config files", "*.conf"), ("All files", "*.*")))
    root.destroy()
    return filename

def loadsc():
# loads the six circle with a configuration file
    config_file = select_file_to_open(initialdir="/IXS2/data/", filetypes=(("Config files", "*.conf"),))
    if config_file:
        # print(f"Loading six-circle configuration from {config_file}")
        hklps.sc.load(config_file)


def savesc():
# saves the six circle configuration to a file
    config_file = select_file_to_save(initialdir="/IXS2/data/")
    if config_file:
        print(f"Saving six-circle configuration to {config_file}")
        hklps.sc.save(config_file)


def sc_init():
# initializes the six circle with the configuration file
    # config_file = select_file(initialdir="/IXS2/data/")
    # hklps.sc.load(config_file)
    hklps.sc.wh_off()
    hklps.sc.setfrozen(345)
    hklps.sc.freeze(0,-0.401,0.401)
    hklps.sc.mv(mu=-0.401, gam=0.401)
    hklps.sc.wh()
#    update_spec_signals()

sc_init()

# def hkl_positions():
# # prints the positions of the pseudomotors (H, K, L) and real motors (Tth, Th, Chi, Phi)
#     sc.wh_off()
# # Print pseudomotor positions (H, K, L)
#     print("Pseudopositions:", hklps.position)

#     # Print real motor positions (Tth, Th, Chi, Phi)
#     print("Real positions:", hklps.real_position)
#     print()
#     # Or, for individual values:
#     print(f"H: {hklps.position.H}")
#     print(f"K: {hklps.position.K}")
#     print(f"L: {hklps.position.L}")
#     print()
#     print(f"Tth: {hklps.real_position.tth}")
#     print(f"Th: {hklps.real_position.th}")
#     print(f"Chi: {hklps.real_position.chi}")
#     print(f"Phi: {hklps.real_position.phi}")

