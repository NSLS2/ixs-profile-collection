'''
from bluesky.callbacks.olog import OlogCallback
from bluesky.callbacks import LivePlot


from suitcase.spec import DocumentToSpec

spec_cb = DocumentToSpec('/tmp/spec1.spec')
RE.subscribe(spec_cb)
'''

import tkinter as tk
from tkinter import filedialog


def wh():
    # updates the H, K, L pseudomotors from the current position of the sample
    hklps.sc.wh()


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


# def select_spec_file(initialdir="."):
#     # Opens a file dialog to select a directory
#     root = tk.Tk()
#     root.withdraw()
#     filename = filedialog.asksaveasfilename(initialdir=initialdir, filetypes=(("SPEC files", "*.spec"), ("All files", "*.*")))
#     root.destroy()
#     return filename


# def newfile():
# # opens a new spec-file
#     config_file = select_spec_file(initialdir="/IXS2/data/")

#     if config_file:
#         # Split path and extension safely
#         directory, prefix = os.path.split(config_file)
#         prefix = os.path.splitext(prefix)[0]  # remove .spec if present

#         # Assign to both RE.md and spec_factory
#         RE.md['spec_file'] = config_file
#         spec_factory.directory = directory
#         spec_factory.prefix = prefix

#         print(f"New SPEC file selected:\n Directory: {directory}\n Prefix: {prefix}")
#     else:
#         print("File selection cancelled.")

def newfile():
    """
    Interactive function to view and modify SPEC-file directory and file name.
   """
    while True:
        print("\n--- Current File Settings ---")
        print(f"Directory: {spec_factory.directory}")
        print(f"File name: {spec_factory.prefix}")
        print("------------------------------")
        print("Options:")
        print("1. Change directory path")
        print("2. Change file name")
        print("3. Exit")
        
        choice = input("Enter your choice (1–3): ").strip()
        
        if choice == "1":
            new_path = input("Enter new directory path: ").strip()
            if new_path:
                spec_factory.directory = new_path
                print(f"Directory updated to: {spec_factory.directory}")
        elif choice == "2":
            new_name = input("Enter new file name: ").strip()
            if new_name:
                spec_factory.prefix = new_name
                print(f"File name updated to: {spec_factory.prefix}")
        elif choice == "3":
            print("Exiting file settings manager.")
            break
        else:
            print("Invalid input. Please enter 1, 2, or 3.")
        
        RE.md['spec_file'] = os.path.join(spec_factory.directory, spec_factory.prefix + ".spec")

    # return filepath, filename


def sc_init():
# initializes the six circle with the configuration file
    # config_file = select_file(initialdir="/IXS2/data/")
    # hklps.sc.load(config_file)
    hklps.sc.__doc__()
    hklps.sc.setprecision(6)
    hklps.sc.wh_off()
    hklps.sc.setfrozen(345)
    hklps.sc.freeze(0,-0.401,0.401)
    hklps.sc.mv(mu=-0.401, gam=0.401)
    hklps.sc.wh()
#    update_spec_signals()

def print_all(H, K, L):
    # prints all angles and positions
    # H, K, L = 1, 1, 1
    flag, pos = ca_s(H,K,L)
    if flag:
        print(f"Number of solutions: {len(pos)}")
        print(' ')
        Mn = ['tth', 'th', 'chi', 'phi', 'mu', 'gam', 'sa', 'omega', 'azimuth', 'alpha', 'beta']
        for n in range(0,len(pos)):
            print(f"Solution {n+1}:")
            res = pos[n]
            caTTH, caTH, caCHI, caPHI, caMU, caGAM, caSA, caOMEGA, caAZIMUTH, caALPHA, caBETA = pos[n]
            caABSQ = hklps.sc.scbasic.Q_length(hklps.sc.LAMBDA,abs(caSA/2))
            # print('')
            print ('H K L =  {0:.{3}f}  {1:.{3}f}  {2:.{3}f}'.format(H,K,L,3))
            print ('|Q| = {0:.3f} nm-1  SA = {1:.{3}f} deg  at  LAMBDA = {2:.{4}f} A'.format(caABSQ*10,caSA,hklps.sc.LAMBDA,3,5))
            print ('AZ = ({0}, {1}, {2})  AZIMUTH = {3:.{6}f} deg  ALPHA = {4:.{6}f}  BETA = {5:.{6}f}'.format(hklps.sc.g_haz,hklps.sc.g_kaz,hklps.sc.g_laz,caAZIMUTH,caALPHA,caBETA,3))
            print ('Omega = th-tth/2 = {0:.{1}f}'.format(caOMEGA,3))
            print ('')
            strfmt = ('{:>'+str(9)+'}')*6
            strprt = ('tth','th','chi','phi','mu','gam')
            print (strfmt.format(*strprt))
            posfmt = ('{:>'+str(9)+'.'+str(4)+'f}')*6
            posprt = (caTTH,caTH,caCHI,caPHI,caMU,caGAM)
            print (posfmt.format(*posprt))
            print ('')


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

