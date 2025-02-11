# Viewer for sixcircle_rqd

import sixcircle
from sixcircle import *
import sixcircle_rqd
from sixcircle_rqd import *

import numpy as np
from numpy.linalg import norm
import math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib.text import Annotation
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.text import Annotation

class Annotation3D(Annotation):
    '''Annotate the point xyz with text s'''

    def __init__(self, s, xyz, *args, **kwargs):
        Annotation.__init__(self,s, xy=(0,0), *args, **kwargs)
        self._verts3d = xyz        

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.xy=(xs,ys)
        Annotation.draw(self, renderer)

def annotate3D(ax, s, *args, **kwargs):
    '''add anotation text s to to Axes3d ax'''
    tag = Annotation3D(s, *args, **kwargs)
    ax.add_artist(tag)        

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


# Calculate H, K, L of every analyzer with viewer
# Usage:  vca6(H,K,L)  or  vca6(H,K,L,Href,Kref,Lref)
# Reference vector (Href, Kref, Lref), usually is a nearby gamma point
sixcircle.wdesc('vca6(H,K,L)','Finds analyzer Q vectors when the center of the arm at (H,K,L)')
sixcircle.wdesc('vca6(H,K,L,Href,Kref,Lref)','Finds analyzer q = (H,K,L)-(Href,Kref,Lref) when the center of the arm at (H,K,L)')
def vca6(*args):
    for value in args:
        if (type(value) in [int,float]) == False:
            print ('\nInvalid argument: {0}\n'.format(value))
    if len(args) == 3:
        caH, caK, caL = args
    elif len(args) == 6:
        caH, caK, caL, Href, Kref, Lref = args
    else:
        print ('\nUsage:  vca6(H,K,L)  or  vca6(H,K,L,Href,Kref,Lref)\n')
        return
    # Calculate positions
    flagca, pos = ca_s(caH,caK,caL)
    if flagca == False:
        print ('Error: Impossible reflection for present frozen!')
        return
    # Record current positions for ease of coming back after calculation
    o_tth, o_th, o_chi, o_phi, o_mu, o_gam = (sixcircle.TTH, sixcircle.TH, sixcircle.CHI, sixcircle.PHI, sixcircle.MU, sixcircle.GAM)
    # Calculate for the first set of pos
    ca_tth, ca_th, ca_chi, ca_phi, ca_mu, ca_gam, ca_sa, ca_omega, ca_azimuth, ca_alpha, ca_beta = pos[0]
    # Record current condition of wh_on (wh_off); turn off wh() after mv()
    o_FLAG_WH = sixcircle.FLAG_WH
    sixcircle.FLAG_WH = False
    mv(tth=ca_tth,th=ca_th,chi=ca_chi,phi=ca_phi,mu=ca_mu,gam=ca_gam)
    A_str = [['0' for xi in range(0,x_n)] for zi in range(0,z_n)]
    A_Q = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    A_q = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    A_ABSQ = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    A_absq = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    A_dQH = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    A_dQV = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    A_dQ = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    Q = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    q = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]    
    for zi in range(0,z_n):
        for xi in range(0,x_n):
            # Give analyzer name
            if flag_bl == 43:
                if zi <= 2:
                    A_str[zi][xi] = 'a' + str(11*zi+6+(xi-x_cen)).zfill(2)
                elif zi >= 3:
                    A_str[zi][xi] = 'a' + str(11*zi+6+(xi-x_cen)-1).zfill(2)
            if flag_bl == 35:
                if zi == 0:
                    A_str[zi][xi] = 'a' + str(5+xi).zfill(2)
                elif zi == 1:
                    A_str[zi][xi] = 'a' + str(1+xi).zfill(2)
                elif zi == 2:               
                    A_str[zi][xi] = 'a' + str(9+xi).zfill(2)
            # Calculate H, K, L at center of each analyzer (A_gam, A_tth)
            A_z = z_off + z_spac * (z_cen - zi)
            A_x = x_off + x_spac * (xi - x_cen)
            A_tth = ca_tth + math.degrees(math.atan2(A_x, a_radi))
            A_gam = ca_gam + math.degrees(math.atan2(A_z, (a_radi**2 + A_x**2)**0.5))
            mv(tth=A_tth, gam=A_gam)
            A_Q[zi][xi] = [sixcircle.H, sixcircle.K, sixcircle.L]
            # Q length, unit A-1
            A_ABSQ[zi][xi] = np.linalg.norm(np.dot(sixcircle.M_B, np.array([[sixcircle.H],[sixcircle.K],[sixcircle.L]])))
            if len(args)==6:
                A_q[zi][xi] = [sixcircle.H-Href, sixcircle.K-Kref, sixcircle.L-Lref]
                # Q length, unit A-1                
                A_absq[zi][xi] = np.linalg.norm(np.dot(sixcircle.M_B, np.array([[sixcircle.H-Href],[sixcircle.K-Kref],[sixcircle.L-Lref]])))
            # Calculate H, K, L at four edges of each analyzer
            A_x_left = A_x + min(agaph*a_radi/sh_radi,ah_size)/2
            A_x_right = A_x - min(agaph*a_radi/sh_radi,ah_size)/2
            A_z_up = A_z + min(agapv*a_radi/sv_radi,av_size)/2
            A_z_low = A_z - min(agapv*a_radi/sv_radi,av_size)/2
            A_tth_left = ca_tth + math.degrees(math.atan2(A_x_left, a_radi))
            A_tth_right = ca_tth + math.degrees(math.atan2(A_x_right, a_radi))
            A_gam_up = ca_gam + math.degrees(math.atan2(A_z_up, (a_radi**2 + A_x**2)**0.5))
            A_gam_low = ca_gam + math.degrees(math.atan2(A_z_low, (a_radi**2 + A_x**2)**0.5))
            # Upper edge
            mv(tth=A_tth, gam=A_gam_up)
            A_Q_up = np.array([sixcircle.H, sixcircle.K, sixcircle.L])
            # Lower edge
            mv(tth=A_tth, gam=A_gam_low)
            A_Q_low = np.array([sixcircle.H, sixcircle.K, sixcircle.L])
            # Left edge
            mv(tth=A_tth_left, gam=A_gam)
            A_Q_left = np.array([sixcircle.H, sixcircle.K, sixcircle.L])
            # Right edge
            mv(tth=A_tth_right, gam=A_gam)
            A_Q_right = np.array([sixcircle.H, sixcircle.K, sixcircle.L])
            # Horizontal dQ
            A_dQH[zi][xi] = A_Q_left - A_Q_right
            # Vertical dQ
            A_dQV[zi][xi] = A_Q_low - A_Q_up
            # dQ
            A_dQ[zi][xi] = (A_dQH[zi][xi]**2 + A_dQV[zi][xi]**2)**0.5
    # Average horizontal dQ and vertical dQ for 28 positions
    aver_dQH = np.array([0,0,0])
    aver_dQV = np.array([0,0,0])
    for zi in range(0,z_n):
        for xi in range(0,x_n):
            aver_dQH = aver_dQH + A_dQH[zi][xi]
            aver_dQV = aver_dQV + A_dQV[zi][xi]
    aver_dQH = aver_dQH / (z_n*x_n)
    aver_dQV = aver_dQV / (z_n*x_n)
    print('')
    print('Q: ({0:.{3}f}  {1:.{3}f}  {2:.{3}f})'.format(caH,caK,caL,sixcircle.PRE), end='')
    print('    at tth={0:.{8}f}, th={1:.{8}f}, chi={2:.{8}f}, phi={3:.{8}f}, mu={4:.{8}f}, gam={5:.{8}f}  H={6:.1f}  V={7:.1f}'.format(ca_tth,ca_th,ca_chi,ca_phi,ca_mu,ca_gam,agaph,agapv,sixcircle.PRE))
    print('')
    latprt = (sixcircle.g_sample,sixcircle.g_aa,sixcircle.g_bb,sixcircle.g_cc,sixcircle.g_al,sixcircle.g_be,sixcircle.g_ga)
    print('Sample {0}    a/b/c {1:.{7}f}/{2:.{7}f}/{3:.{7}f}    alpha/beta/gamma {4:.{7}f}/{5:.{7}f}/{6:.{7}f}'.format(*latprt,sixcircle.PRE))
    infprt = (sixcircle.LAMBDA,sixcircle.g_frozen,sixcircle.g_haz,sixcircle.g_kaz,sixcircle.g_laz,ca_alpha,ca_beta)
    print('Wavelength {0:.{7}f}    frozen={1}    AZ ({2:.{8}f}, {3:.{8}f}, {4:.{8}f})    ALPHA={5:.{8}f}  BETA={6:.{8}f}'.format(*infprt,sixcircle.PRE+2,sixcircle.PRE))
    print('Or0: ({0:.{3}f}, {1:.{3}f}, {2:.{3}f})'.format(sixcircle.g_h0,sixcircle.g_k0,sixcircle.g_l0,sixcircle.PRE), end='')
    or0prt = (sixcircle.g_u00,sixcircle.g_u01,sixcircle.g_u02,sixcircle.g_u03,sixcircle.g_u04,sixcircle.g_u05)
    print('    at tth={0:.{6}f}, th={1:.{6}f}, chi={2:.{6}f}, phi={3:.{6}f}, mu={4:.{6}f}, gam={5:.{6}f}  '.format(*or0prt,sixcircle.PRE))
    print('Or1: ({0:.{3}f}, {1:.{3}f}, {2:.{3}f})'.format(sixcircle.g_h1,sixcircle.g_k1,sixcircle.g_l1,sixcircle.PRE), end='')
    or1prt = (sixcircle.g_u10,sixcircle.g_u11,sixcircle.g_u12,sixcircle.g_u13,sixcircle.g_u14,sixcircle.g_u15)
    print('    at tth={0:.{6}f}, th={1:.{6}f}, chi={2:.{6}f}, phi={3:.{6}f}, mu={4:.{6}f}, gam={5:.{6}f}'.format(*or1prt,sixcircle.PRE))
    print('')
    if len(args) == 6: print('         h         k         l')
    if len(args) == 3: print('         H         K         L')    
    for zi in range(0,z_n):
        for xi in range(0,x_n):
            if len(args) == 6:
                afmt = '{0}: ({1:>{11}.{12}f}, {2:>{11}.{12}f}, {3:>{11}.{12}f})  |q|={4:>7.3f} nm-1  dq:({5:>{11}.{12}f}, {6:>{11}.{12}f}, {7:>{11}.{12}f})  Qtot:({8:>{11}.{12}f}, {9:>{11}.{12}f}, {10:>{11}.{12}f})'
                aprt = (A_str[zi][xi],A_q[zi][xi][0],A_q[zi][xi][1],A_q[zi][xi][2],A_absq[zi][xi]*10,A_dQ[zi][xi][0],A_dQ[zi][xi][1],A_dQ[zi][xi][2],A_Q[zi][xi][0],A_Q[zi][xi][1],A_Q[zi][xi][2])
                print(afmt.format(*aprt,sixcircle.PRE+4,sixcircle.PRE))
            if len(args) == 3:
                afmt = '{0}: ({1:>{8}.{9}f}, {2:>{8}.{9}f}, {3:>{8}.{9}f})  |Q|={4:>7.3f} nm-1  dq:({5:>{8}.{9}f}, {6:>{8}.{9}f}, {7:>{8}.{9}f})'
                aprt = (A_str[zi][xi],A_Q[zi][xi][0],A_Q[zi][xi][1],A_Q[zi][xi][2],A_ABSQ[zi][xi]*10,A_dQ[zi][xi][0],A_dQ[zi][xi][1],A_dQ[zi][xi][2],)
                print(afmt.format(*aprt,sixcircle.PRE+4,sixcircle.PRE))
        print('')
    dqprt = (agaph,aver_dQH[0],aver_dQH[1],aver_dQH[2],agapv,aver_dQV[0],aver_dQV[1],aver_dQV[2])
    print('Av.  dq  H({0:.1f}): ({1:.{8}f},{2:.{8}f},{3:.{8}f})  V({4:.1f}): ({5:.{8}f},{6:.{8}f},{7:.{8}f})'.format(*dqprt,sixcircle.PRE))
    # Come back to original positions
    mv(tth=o_tth, th=o_th, chi=o_chi, phi=o_phi, mu=o_mu, gam=o_gam)
    # Come back to original wh_on (wh_off)
    sixcircle.FLAG_WH = o_FLAG_WH
    # Create document gpi.hkl_pos
    print ('')
    print ('HKL values to:  gpi.hkl_pos')
    try:
        with open('gpi.hkl_pos', 'w') as f:
            for zi in range(0,z_n):
                for xi in range(0,x_n):
                    f.write('hkl_{0} = "({1:.{4}f}, {2:.{4}f}, {3:.{4}f})"\n'.format(A_str[zi][xi],A_Q[zi][xi][0],A_Q[zi][xi][1],A_Q[zi][xi][2],sixcircle.PRE))        
            f.write('\n')
            f.write('dq_Hav({0:.1f}) = "({1:.{4}f}, {2:.{4}f}, {3:.{4}f})"\n'.format(agaph,aver_dQH[0],aver_dQH[1],aver_dQH[2],sixcircle.PRE))
            f.write('dq_Vav({0:.1f}) = "({1:.{4}f}, {2:.{4}f}, {3:.{4}f})"'.format(agapv,aver_dQV[0],aver_dQV[1],aver_dQV[2],sixcircle.PRE))
    except:
        print ('\nError in writing gpi.hkl_pos')
    print ('')
    print ('Command(BL43LXU):      ', end='')    
    print ('mv tth {0:.{8}f} th {1:.{8}f} chi {2:.{8}f} phi {3:.{8}f} agaph {6:.1f} agapv {7:.1f}'.format(ca_tth,ca_th,ca_chi,ca_phi,ca_mu,ca_gam,agaph,agapv,sixcircle.PRE))
    print ()
    if len(args) == 6: print('         qx        qy        qz')
    if len(args) == 3: print('         Qx        Qy        Qz')    
    for zi in range(0,z_n):
        for xi in range(0,x_n):
            if len(args) == 6:
                afmt = '{0}: ({1:>{11}.{12}f}, {2:>{11}.{12}f}, {3:>{11}.{12}f})  |q|={4:>7.3f} nm-1  dq:({5:>{11}.{12}f}, {6:>{11}.{12}f}, {7:>{11}.{12}f})  Qtot:({8:>{11}.{12}f}, {9:>{11}.{12}f}, {10:>{11}.{12}f})'
                q[zi][xi] = UBT(A_q[zi][xi][0],A_q[zi][xi][1],A_q[zi][xi][2])
                aprt = (A_str[zi][xi],q[zi][xi][0],q[zi][xi][1],q[zi][xi][2],A_absq[zi][xi]*10,A_dQ[zi][xi][0],A_dQ[zi][xi][1],A_dQ[zi][xi][2],A_Q[zi][xi][0],A_Q[zi][xi][1],A_Q[zi][xi][2])
                print(afmt.format(*aprt,sixcircle.PRE+4,sixcircle.PRE))
            if len(args) == 3:
                afmt = '{0}: ({1:>{8}.{9}f}, {2:>{8}.{9}f}, {3:>{8}.{9}f})  |Q|={4:>7.3f} nm-1  dq:({5:>{8}.{9}f}, {6:>{8}.{9}f}, {7:>{8}.{9}f})'
                Q[zi][xi] = UBT(A_Q[zi][xi][0],A_Q[zi][xi][1],A_Q[zi][xi][2])
                aprt = (A_str[zi][xi],Q[zi][xi][0],Q[zi][xi][1],Q[zi][xi][2],A_ABSQ[zi][xi]*10,A_dQ[zi][xi][0],A_dQ[zi][xi][1],A_dQ[zi][xi][2],)
                print(afmt.format(*aprt,sixcircle.PRE+4,sixcircle.PRE))
        print('')
    # 
    if len(args) == 6:
        plt_qxyz(q,args)        
        #plt_hkl_cubic(A_q,args)
    if len(args) == 3:
        plt_qxyz(Q,args)
        #plt_hkl_cubic(A_Q,args)
    
    
def plt_qxyz(q,args):
    X0 = [q[0][xi][0] for xi in range(0,x_n)]; Y0 = [q[0][xi][1] for xi in range(0,x_n)]; Z0 = [q[0][xi][2] for xi in range(0,x_n)]   
    X1 = [q[1][xi][0] for xi in range(0,x_n)]; Y1 = [q[1][xi][1] for xi in range(0,x_n)]; Z1 = [q[1][xi][2] for xi in range(0,x_n)]   
    X2 = [q[2][xi][0] for xi in range(0,x_n)]; Y2 = [q[2][xi][1] for xi in range(0,x_n)]; Z2 = [q[2][xi][2] for xi in range(0,x_n)]   
    X3 = [q[3][xi][0] for xi in range(0,x_n)]; Y3 = [q[3][xi][1] for xi in range(0,x_n)]; Z3 = [q[3][xi][2] for xi in range(0,x_n)]   
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d', aspect='auto')
    ax.axis("off")
    ax.scatter(X0, Y0, Z0); ax.scatter(X1, Y1, Z1); ax.scatter(X2, Y2, Z2); ax.scatter(X3, Y3, Z3)    
    if len(args) == 6: ax.set_xlabel('$q_x$', fontsize=15);ax.set_ylabel('$q_y$', fontsize=15);ax.set_zlabel('$q_z$', fontsize=15)
    if len(args) == 3: ax.set_xlabel('$Q_x$', fontsize=15);ax.set_ylabel('$Q_y$', fontsize=15);ax.set_zlabel('$Q_z$', fontsize=15)        
    latprt = (sixcircle.g_sample,sixcircle.g_aa,sixcircle.g_bb,sixcircle.g_cc,sixcircle.g_al,sixcircle.g_be,sixcircle.g_ga)
    ax.set_title('Sample {0}    a/b/c {1:.{7}f}/{2:.{7}f}/{3:.{7}f}    alpha/beta/gamma {4:.{7}f}/{5:.{7}f}/{6:.{7}f}'.format(*latprt,sixcircle.PRE),fontsize=8)
    q100 = UBT(1,0,0)
    q010 = UBT(0,1,0)
    q001 = UBT(0,0,1)
    q100 /= norm(q100); q010 /= norm(q010); q001 /= norm(q001)
    q100 *= 3.;q010 *= 3.;q001 *= 3;
    O = np.array([0, 0, 0])
    a = Arrow3D([0,q100[0]], [0,q100[1]], [0,q100[2]], mutation_scale=1, lw=1, arrowstyle="-|>", color="b")
    b = Arrow3D([0,q010[0]], [0,q010[1]], [0,q010[2]], mutation_scale=1, lw=1, arrowstyle="-|>", color="y")
    c = Arrow3D([0,q001[0]], [0,q001[1]], [0,q001[2]], mutation_scale=1, lw=1, arrowstyle="-|>", color="g")
    ax.add_artist(a); ax.add_artist(b); ax.add_artist(c)
    annotate3D(ax, s=str('[100]'), xyz=q100, fontsize=7, xytext=(-3,3), textcoords='offset points', ha='center',va='center')
    annotate3D(ax, s=str('[010]'), xyz=q010, fontsize=7, xytext=(-3,3), textcoords='offset points', ha='center',va='center')
    annotate3D(ax, s=str('[001]'), xyz=q001, fontsize=7, xytext=(-3,3), textcoords='offset points', ha='center',va='center')  
    plt.show() 

def plt_hkl_cubic(h,args):
    X0 = [h[0][xi][0] for xi in range(0,x_n)]; Y0 = [h[0][xi][1] for xi in range(0,x_n)]; Z0 = [h[0][xi][2] for xi in range(0,x_n)]   
    X1 = [h[1][xi][0] for xi in range(0,x_n)]; Y1 = [h[1][xi][1] for xi in range(0,x_n)]; Z1 = [h[1][xi][2] for xi in range(0,x_n)]   
    X2 = [h[2][xi][0] for xi in range(0,x_n)]; Y2 = [h[2][xi][1] for xi in range(0,x_n)]; Z2 = [h[2][xi][2] for xi in range(0,x_n)]   
    X3 = [h[3][xi][0] for xi in range(0,x_n)]; Y3 = [h[3][xi][1] for xi in range(0,x_n)]; Z3 = [h[3][xi][2] for xi in range(0,x_n)]   
    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d', aspect='auto')
    ax.scatter(X0, Y0, Z0); ax.scatter(X1, Y1, Z1); ax.scatter(X2, Y2, Z2); ax.scatter(X3, Y3, Z3)    
    if len(args) == 6: ax.set_xlabel('$h$', fontsize=15);ax.set_ylabel('$k$', fontsize=15);ax.set_zlabel('$l$', fontsize=15)
    if len(args) == 3: ax.set_xlabel('$H$', fontsize=15);ax.set_ylabel('$K$', fontsize=15);ax.set_zlabel('$L$', fontsize=15)        
    latprt = (sixcircle.g_sample,sixcircle.g_aa,sixcircle.g_bb,sixcircle.g_cc,sixcircle.g_al,sixcircle.g_be,sixcircle.g_ga)
    ax.set_title('Sample {0}    a/b/c {1:.{7}f}/{2:.{7}f}/{3:.{7}f}    alpha/beta/gamma {4:.{7}f}/{5:.{7}f}/{6:.{7}f}'.format(*latprt,sixcircle.PRE),fontsize=8)
    plt.show()
    #plt.show(block=False)     


# Transform reciprocal lattice unit to laboratory coordinate by UB matrics
# output unit: nm^-1
def UBT(h,k,l):
    q = sixcircle.M_U@sixcircle.M_B@[h,k,l]; q = 10.*q
    return q
    
