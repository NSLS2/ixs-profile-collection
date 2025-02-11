# sixcircle_rqd_view: analyzer viewer for diffractometers at BL35XU, BL43LXU, SPring-8
# Required additional documents (for included codess): 0- sixcircle_rqd.py 1- sixcircle.py, 2- scbasic.py, 3- ini.conf, 4- BL43XU_CONST.mac, 5- BL35XU_CONST.mac

# Developed with Python 3.7.3, numpy 1.16.4, scipy 1.3.0
# Materials Dynamics Laboratory, RIKEN SPring-8 Center
# Ver 1.52, December 2020

# Authors: Daisuke Ishikawa & Alfred Q. R. Baron
# Contact: disikawa@spring8.or.jp, baron@spring8.or.jp

# If this code is used independently of work at BL43, please reference the writeup ("Open-source Python software for six-circle diffraction with an inelastic x-ray scattering (IXS) spectrometer" by Wenyang ZHAO and Alfred Q.R. BARON, unpublished, available at https://beamline.harima.riken.jp/bl43lxu)

# In typical use related to experimental work at BL43LXU, it is enough to reference a BL publication.  However, if this code is used extensively, then please include a separate reference as above.

# Viewer for sixcircle_rqd

import sixcircle
from sixcircle import *
import sixcircle_rqd
from sixcircle_rqd import *

import numpy as np
import math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib.text import Annotation
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.text import Annotation

print ('Loading sixcircle_rqd_view.py - for IXS users at SPring-8 - D. Ishikawa and A. Baron')

#########################################################################################################################################
# Reciporcal space definitions - used in both 2D and 3D plotting
#########################################################################################################################################
def cross3x3(v1,v2):
#manual cross product because numpy is SLOW for 3-vectors
    r0 = v1[1]*v2[2] - v1[2]*v2[1]
    r1 = v1[2]*v2[0] - v1[0]*v2[2]
    r2 = v1[0]*v2[1] - v1[1]*v2[0]
    cp = [r0,r1,r2]
    return cp

# find catesian components of real space and recprocal space vectors using the orientation matrix (a,b,c,alpha,beta,gamma)
def get_cartesian_vectors():
    cosalp = np.cos(np.radians(g_al)) ; sinalp = np.sin(np.radians(g_al))
    cosbet = np.cos(np.radians(g_be)) ; sinbet = np.sin(np.radians(g_be))
    cosgam = np.cos(np.radians(g_ga)) ; singam = np.sin(np.radians(g_ga))
    av = [g_aa, 0. , 0.]
    bv = [g_bb * cosgam , g_bb * singam , 0.]
    cv = [g_cc * cosbet , g_cc*(cosalp-cosbet*cosgam)/singam, g_cc * np.sqrt(1. - cosbet**2 - ((cosalp-cosbet*cosgam)/singam)**2)]
    astar = cross3x3(bv,cv)
    bstar = cross3x3(cv,av)
    cstar = cross3x3(av,bv)
    vol = g_aa*g_bb*g_cc * np.sqrt(1+2.*cosalp*cosbet*cosgam-cosalp**2-cosbet**2-cosgam**2)
    sf = 2.*np.pi/vol
    astar = [sf*astar[0],sf*astar[1],sf*astar[2]]; bstar = [sf*bstar[0],sf*bstar[1],sf*bstar[2]]; cstar = [sf*cstar[0],sf*cstar[1],sf*cstar[2]]
    #print('         a=',av,'        b=',bv,'        c=',cv)
    #print('     astar=',astar,'    bstar=',bstar,'    cstar=',cstar)
    return av,bv,cv,astar,bstar,cstar

def hkl_to_cartesian(vhkl,Hcart,Kcart,Lcart):
    vcart = [0.,0.,0.]
    for i in range(3):
        vcart[i] = vhkl[0]*Hcart[i] + vhkl[1]*Kcart[i] + vhkl[2]*Lcart[i]
    return vcart

def get_analyzer_ranges(Q,dQH,dQV):
    qmin = [0,0,0]; qmax= [0,0,0]
    for j in range(3):
        qmin[j] = Q[0][0][j]; qmax[j] = Q[0][0][j]
    for j in range(3):
        for iz in range(0,z_n):
            for ix in range(0,x_n):
                c = Q[iz][ix][j]; d1 = dQH[iz][ix][j]; d2 = dQV[iz][ix][j]
                ll = [c+d1+d2,c+d1-d2,c-d1+d2,c-d1-d2,qmin[j],qmax[j]]
                qmin[j] = min(ll)
                qmax[j] = max(ll)
    return qmin,qmax

#########################################################################################################################################
# 3D plotting definitions
#########################################################################################################################################
sixcircle.wdesc('ca6_view3D(H,K,L)','3D plot of analyzer Qs')
def ca6_view3D(*args):
    ca6(*args)
    # transfer variables
    A_str = sixcircle_rqd.A_str;  A_Q = sixcircle_rqd.A_Q ; A_q = sixcircle_rqd.A_q
    A_ABSQ=sixcircle_rqd.A_ABSQ; A_absq = sixcircle_rqd.A_absq
    A_dQH=sixcircle_rqd.A_dQH; A_dQV=sixcircle_rqd.A_dQV; A_dQ=sixcircle_rqd.A_dQ
    global x_n,z_n
    x_n = sixcircle_rqd.x_n; z_n = sixcircle_rqd.z_n
    # convert to orthonormal cartesian coordinates
    (av,bv,cv,Hcart,Kcart,Lcart) = get_cartesian_vectors()
    Qcart = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    qcart = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    dQcartH = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    dQcartV = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    for zi in range(0,z_n):
        for xi in range(0,x_n):
            if len(args) == 6:
                qcart[zi][xi] = hkl_to_cartesian(A_q[zi][xi],Hcart,Kcart,Lcart)
            Qcart[zi][xi] = hkl_to_cartesian(A_Q[zi][xi],Hcart,Kcart,Lcart)
            dQcartV[zi][xi] = hkl_to_cartesian(A_dQV[zi][xi],Hcart,Kcart,Lcart)
            dQcartH[zi][xi] = hkl_to_cartesian(A_dQH[zi][xi],Hcart,Kcart,Lcart)
    print("\nMaking 3D plot.")
    make_analyzer_figure_3D(Hcart,Kcart,Lcart,Qcart,dQcartH,dQcartV,A_str)

def make_analyzer_figure_3D(H,K,L,Q,dQH,dQV,A_str):
    global ax
    #print("H = ",H,"    K=",K,"    L=",L)
    fig = plt.figure(figsize=[6,6],linewidth=3.,frameon=None)
    ax = fig.gca(projection='3d')
    qmin=[0,0,0]; qmax=[0,0,0]
    qmin,qmax = get_analyzer_ranges(Q,dQH,dQV)
    ax.set_xlim(qmin[0],qmax[0])
    ax.set_ylim(qmin[1],qmax[1])
    ax.set_zlim(qmin[2],qmax[2])
    ax.set_axis_off()
    # get arm center
    arm_center=[0,0,0]
    for zi in range(0,z_n):
        for xi in range(0,x_n):
            if sixcircle_rqd.flag_bl == 43 :
                if A_str[zi][xi] == 'a06' : arm_center = Q[zi][xi]
            if sixcircle_rqd.flag_bl == 35 :
                if A_str[zi][xi] == 'a02' or A_str[zi][xi]=='a03' :
                    for j in range(3) : arm_center[j] = arm_center[j]+Q[zi][xi]/2.
    #print("arm_center = ",arm_center)
    plot_directions_3D(arm_center,[H,K,L])

    for zi in range(0,z_n):
        for xi in range(0,x_n):
            #print(A_str[zi][xi],Q[zi][xi])
            plot_one_analyzer_3D(Q[zi][xi],dQH[zi][xi],dQV[zi][xi],A_str[zi][xi])
    plot_grid_3D(H,K,L,[qmin[0]-1,qmax[0]+1],[qmin[1]-1,qmax[1]+1],[qmin[2]-1,qmax[2]+1])
    plot_latticepoints_3D(H,K,L,[qmin[0]-1,qmax[0]+1],[qmin[1]-1,qmax[1]+1],[qmin[2]-1,qmax[2]+1])
    plt.show()

def plot_one_analyzer_3D(center,wid1,wid2,alabel):
    hwid1 = [i/2. for i in wid1] ;    hwid2 = [i/2. for i in wid2]  # half-width
    xlist = 5*[center[0]] ;     ylist = 5*[center[1]];     zlist = 5*[center[2]]  # define the arrays
    xlist[0] = center[0] - hwid1[0] - hwid2[0];   ylist[0] = center[1] - hwid1[1] - hwid2[1];   zlist[0] = center[2] - hwid1[2] - hwid2[2]
    xlist[1] = center[0] - hwid1[0] + hwid2[0];   ylist[1] = center[1] - hwid1[1] + hwid2[1];   zlist[1] = center[2] - hwid1[2] + hwid2[2]
    xlist[2] = center[0] + hwid1[0] + hwid2[0];   ylist[2] = center[1] + hwid1[1] + hwid2[1];   zlist[2] = center[2] + hwid1[2] + hwid2[2]
    xlist[3] = center[0] + hwid1[0] - hwid2[0];   ylist[3] = center[1] + hwid1[1] - hwid2[1];   zlist[3] = center[2] + hwid1[2] - hwid2[2]
    xlist[4] = center[0] - hwid1[0] - hwid2[0];   ylist[4] = center[1] - hwid1[1] - hwid2[1];   zlist[4] = center[2] - hwid1[2] - hwid2[2]
    ax.plot(xlist,ylist,zlist,'k')  # draw rectangle border for analyzer
    #ax.text(center[0],center[1],center[2],alabel) # label it
    xlist = [0.,center[0]]; ylist = [0.,center[1]]; zlist=[0.,center[2]]
    ax.plot(xlist,ylist,zlist,'k',lw=0.2)
    plt.pause(0.001)

def plot_directions_3D(center,v):
    llw = 3.
    llc = ['r','b','g']
    textstr = ['a*','b*','c*']
    for iv in range(3) :
        endpt = [(center[i]+v[iv][i]) for i in range(3)]
        xlist = [center[0],endpt[0]] ;    ylist = [center[1],endpt[1]] ;   zlist = [center[2],endpt[2]]
        ax.plot(xlist,ylist,zlist,llc[iv],lw=llw)
        ax.text(endpt[0],endpt[1],endpt[2],textstr[iv])
        halfpt = [(center[i]+v[iv][i]/2.) for i in range(3)]
        ax.text(halfpt[0],halfpt[1],halfpt[2],textstr[iv])

def plot_grid_3D(v1,v2,v3,r1,r2,r3):
    llw = 0.3; lls = '-'; lc='b'
    for i in range( int(r1[0]) , int(r1[1]) ) :
        for j in range( int(r2[0]) , int(r2[1]) ) :
            refpt = [ i*v1[0]+j*v2[0], i*v1[1]+j*v2[1], i*v1[2]+j*v2[2] ]
            xlist = [ refpt[0] + r3[0] * v3[0] , refpt[0] + r3[1] *v3[0] ]
            ylist = [ refpt[1] + r3[0] * v3[1] , refpt[1] + r3[1] *v3[1] ]
            zlist = [ refpt[2] + r3[0] * v3[2] , refpt[2] + r3[1] *v3[2] ]
            ax.plot(xlist,ylist,zlist,lc,lw=llw,ls=lls)
    for i in range( int(r2[0]) , int(r2[1]) ) :
        for j in range( int(r3[0]) , int(r3[1]) ) :
            refpt = [ i*v2[0]+j*v3[0], i*v2[1]+j*v3[1], i*v2[2]+j*v3[2] ]
            xlist = [ refpt[0] + r1[0] * v1[0] , refpt[0] + r1[1] *v1[0] ]
            ylist = [ refpt[1] + r1[0] * v1[1] , refpt[1] + r1[1] *v1[1] ]
            zlist = [ refpt[2] + r1[0] * v1[2] , refpt[2] + r1[1] *v1[2] ]
            ax.plot(xlist,ylist,zlist,lc,lw=llw,ls=lls)
    for i in range( int(r1[0]) , int(r1[1]) ) :
        for j in range( int(r3[0]) , int(r3[1]) ) :
            refpt = [ i*v1[0]+j*v3[0], i*v1[1]+j*v3[1], i*v1[2]+j*v3[2] ]
            xlist = [ refpt[0] + r2[0] * v2[0] , refpt[0] + r2[1] *v2[0] ]
            ylist = [ refpt[1] + r2[0] * v2[1] , refpt[1] + r2[1] *v2[1] ]
            zlist = [ refpt[2] + r2[0] * v2[2] , refpt[2] + r2[1] *v2[2] ]
            ax.plot(xlist,ylist,zlist,lc,lw=llw,ls=lls)

def plot_latticepoints_3D(v1,v2,v3,r1,r2,r3):
    llw = 0.3; lls = 'o'; lc='k'
    for i in range( int(r1[0]) , int(r1[1]) ) :
        for j in range( int(r2[0]) , int(r2[1]) ) :
            for k in range( int(r3[0]) , int(r3[1]) ) :
                pt = [ i*v1[0]+j*v2[0]+k*v3[0], i*v1[1]+j*v2[1]+k*v3[1], i*v1[2]+j*v2[2]+k*v3[2] ]
                ax.scatter(pt[0],pt[1],zs=pt[2],s=80,c='k',depthshade=False)
                pstr = "(%d,%d,%d)" %(i,j,k)
                ax.text(pt[0],pt[1],pt[2],pstr)

#########################################################################################################################################
# 2D plotting definitions
#########################################################################################################################################
sixcircle.wdesc('ca6_view2D_HKL(H,K,L)','2D plot of analyzers in HKL')
def ca6_view2D_HKL(*args):
    ca6(*args)
    # transfer variables
    A_str = sixcircle_rqd.A_str;  A_Q = sixcircle_rqd.A_Q ; A_q = sixcircle_rqd.A_q
    A_ABSQ=sixcircle_rqd.A_ABSQ; A_absq = sixcircle_rqd.A_absq
    A_dQH=sixcircle_rqd.A_dQH; A_dQV=sixcircle_rqd.A_dQV; A_dQ=sixcircle_rqd.A_dQ
    global x_n,z_n
    x_n = sixcircle_rqd.x_n; z_n = sixcircle_rqd.z_n
    # convert to orthonormal cartesian coordinates
    print("\nMaking 2D plot vs H K L.")
    make_analyzer_figure_2D_HKL(A_Q,A_q,A_dQH,A_dQV,A_str)

def make_analyzer_figure_2D_HKL(Q,q,dQH,dQV,A_str):
    global plt
    qmin=[0,0,0]; qmax=[0,0,0]
    (qmin,qmax) = get_analyzer_ranges(Q,dQH,dQV)
    pl=['H','K','L']
    fig = plt.figure(figsize=[10,10],linewidth=3.,frameon=None)
    integerlines = True
    fig.add_subplot(223); i=0; j=1;  make_2Dplot(i,j,qmin,qmax,pl,Q,q,dQH,dQV,A_str,integerlines)
    fig.add_subplot(221); i=0; j=2;  make_2Dplot(i,j,qmin,qmax,pl,Q,q,dQH,dQV,A_str,integerlines)
    fig.add_subplot(224); i=2; j=1;  make_2Dplot(i,j,qmin,qmax,pl,Q,q,dQH,dQV,A_str,integerlines)

sixcircle.wdesc('ca6_view2D_xyz(H,K,L)','2D plot of analyzers in orthogonal coords')
def ca6_view2D_xyz(*args):
    ca6(*args)
    # transfer variables
    A_str = sixcircle_rqd.A_str;  A_Q = sixcircle_rqd.A_Q ; A_q = sixcircle_rqd.A_q
    A_ABSQ=sixcircle_rqd.A_ABSQ; A_absq = sixcircle_rqd.A_absq
    A_dQH=sixcircle_rqd.A_dQH; A_dQV=sixcircle_rqd.A_dQV; A_dQ=sixcircle_rqd.A_dQ
    global x_n,z_n
    x_n = sixcircle_rqd.x_n; z_n = sixcircle_rqd.z_n
    # convert to orthonormal cartesian coordinates
    (av,bv,cv,Hcart,Kcart,Lcart) = get_cartesian_vectors()
    Qcart = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    qcart = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    dQcartH = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    dQcartV = [[0 for xi in range(0,x_n)] for zi in range(0,z_n)]
    for zi in range(0,z_n):
        for xi in range(0,x_n):
            if len(args) == 6:
                qcart[zi][xi] = hkl_to_cartesian(A_q[zi][xi],Hcart,Kcart,Lcart)
            Qcart[zi][xi] = hkl_to_cartesian(A_Q[zi][xi],Hcart,Kcart,Lcart)
            dQcartV[zi][xi] = hkl_to_cartesian(A_dQV[zi][xi],Hcart,Kcart,Lcart)
            dQcartH[zi][xi] = hkl_to_cartesian(A_dQH[zi][xi],Hcart,Kcart,Lcart)
    print("\nMaking 2D plot vs Qx,Qy,Qz.")
    make_analyzer_figure_2D_xyz(Hcart,Kcart,Lcart,Qcart,qcart,dQcartH,dQcartV,A_str)

def make_analyzer_figure_2D_xyz(H,K,L,Q,q,dQH,dQV,A_str):
    global plt
    qmin=[0,0,0]; qmax=[0,0,0]
    (qmin,qmax) = get_analyzer_ranges(Q,dQH,dQV)
    pl=['Qx','Qy','Qz']
    fig = plt.figure(figsize=[10,10],linewidth=3.,frameon=None)
    integerlines = False
    fig.add_subplot(223); i=0; j=1;  make_2Dplot(i,j,qmin,qmax,pl,Q,q,dQH,dQV,A_str,integerlines)
    fig.add_subplot(221); i=0; j=2;  make_2Dplot(i,j,qmin,qmax,pl,Q,q,dQH,dQV,A_str,integerlines)
    fig.add_subplot(224); i=2; j=1;  make_2Dplot(i,j,qmin,qmax,pl,Q,q,dQH,dQV,A_str,integerlines)

def make_2Dplot(i,j,qmin,qmax,pl,Q,q,dQH,dQV,A_str,integerlines):
    plt.xlim(qmin[i],qmax[i])
    plt.ylim(qmin[j],qmax[j])
    plt.xlabel(pl[i])
    plt.ylabel(pl[j])
    plt.grid()
    if integerlines :
        locs, labels = plt.xticks()
        for il in range(len(locs)) :
            t = locs[il]
            if t==int(t) :
                xlist = [t,t]; ylist = [qmin[j],qmax[j]]
                plt.plot(xlist,ylist,'k')
        locs, labels = plt.yticks()
        for il in range(len(locs)) :
            t = locs[il]
            if t==int(t) :
                xlist = [qmin[i],qmax[i]]; ylist = [t,t]
                plt.plot(xlist,ylist,'k')
    for zi in range(0,z_n):
        for xi in range(0,x_n):
            plot_one_analyzer_2D(i,j,Q[zi][xi],dQH[zi][xi],dQV[zi][xi],A_str[zi][xi])

def plot_one_analyzer_2D(i1,i2,center,wid1,wid2,alabel):
    hwid1 = [i/2. for i in wid1] ;    hwid2 = [i/2. for i in wid2]  # half-width
    xlist = 5*[center[0]] ;     ylist = 5*[center[1]];     zlist = 5*[center[2]]  # define the arrays
    xlist[0] = center[i1] - hwid1[i1] - hwid2[i1];   ylist[0] = center[i2] - hwid1[i2] - hwid2[i2]
    xlist[1] = center[i1] - hwid1[i1] + hwid2[i1];   ylist[1] = center[i2] - hwid1[i2] + hwid2[i2]
    xlist[2] = center[i1] + hwid1[i1] + hwid2[i1];   ylist[2] = center[i2] + hwid1[i2] + hwid2[i2]
    xlist[3] = center[i1] + hwid1[i1] - hwid2[i1];   ylist[3] = center[i2] + hwid1[i2] - hwid2[i2]
    xlist[4] = center[i1] - hwid1[i1] - hwid2[i1];   ylist[4] = center[i2] - hwid1[i2] - hwid2[i2]
    plt.plot(xlist,ylist,'k')  # draw rectangle border for analyzer
    style = dict(size=8, color='red',va='center',ha='center')
    plt.text(center[i1],center[i2],alabel,**style) # label it
    plt.pause(0.001)





















