#!/usr/bin/env python3
# coding=UTF-8

import sys, time
import os.path
import argparse
import numpy as np
from TimeDomainField import *
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the input file')
parser.add_argument('-ix', help="index x for slicing (default -iy=Ny/2)", dest="ix", type=int)
parser.add_argument('-iy', help="index y for slicing", dest="iy", type=int)
parser.add_argument('-Emax', help="plot range for the E field [V/m]", dest="emax", type=float)
parser.add_argument('-Bmax', help="plot range for the B field [T]", dest="bmax", type=float)

print()
args = parser.parse_args()

infile = args.file
fileOK = os.path.isfile(infile)
if not fileOK:
  print("file not found")
  sys.exit()

Source = TimeDomainField()
# Open the file for reading
print("reading %s" % infile)
Source.read(infile)

print("Nx=%d Ny=%d" % (Source.Nx,Source.Ny))
print("dt=%g NOTS=%d" % (Source.dt, Source.Nt))

if args.ix != None:
    ix = args.ix
    print("slicing ix=%d" % ix)
    A = Source.A[ix,:]
    t0 = Source.t0[ix,:]
    pos = Source.pos[ix,:]
    dx = Source.dyVec
else:
    if args.iy != None:
        iy = args.iy
    else:
        iy = Source.Ny // 2
    print("slicing iy=%d" % iy)
    A = Source.A[:,iy]
    t0 = Source.t0[:,iy]
    pos = Source.pos[:,iy]
    dx = Source.dxVec
print(A.shape)

# get the field components
Ex = A[:,:,0]
Ey = A[:,:,1]
Ez = A[:,:,2]
Bx = A[:,:,3]
By = A[:,:,4]
Bz = A[:,:,5]

# determine the plot coordinates
dt = Source.dt
Nt = Source.Nt
Nx = Ex.shape[1]
t = np.array([range(Nt)*dt+t0S for t0S in t0])
nx = dx / np.linalg.norm(dx)
x0 = np.array([np.dot(p,nx) for p in pos])
x = np.transpose(np.array([x0 for t0S in range(Nt)]))

# determine the plot z-axis levels
if args.emax == None:
    exmax = np.max(np.abs(Ex))
    eymax = np.max(np.abs(Ey))
    ezmax = np.max(np.abs(Ez))
    emax = np.max(np.array([exmax,eymax,ezmax]))
else:
    emax = args.emax
elevels = np.linspace(-emax, emax, num=15)
if args.bmax == None:
    bxmax = np.max(np.abs(Bx))
    bymax = np.max(np.abs(By))
    bzmax = np.max(np.abs(Bz))
    bmax = np.max(np.array([bxmax,bymax,bzmax]))
else:
    bmax = args.bmax
blevels = np.linspace(-bmax, bmax, num=15)

# create the plots
fig, axs = plt.subplots(2, 3, figsize=(12,7),dpi=100,constrained_layout=True)
ax1 = axs[0,0]
ax2 = axs[0,1]
ax3 = axs[0,2]
ax4 = axs[1,0]
ax5 = axs[1,1]
ax6 = axs[1,2]
ax1.contourf(t,x,Ex,elevels,cmap='jet')
ax2.contourf(t,x,Ey,elevels,cmap='jet')
p3 = ax3.contourf(t,x,Ez,elevels,cmap='jet')
ax4.contourf(t,x,Bx,blevels,cmap='jet')
ax5.contourf(t,x,By,blevels,cmap='jet')
p6 = ax6.contourf(t,x,Bz,blevels,cmap='jet')
# annotations
ax1.text(0.85, 0.9, r'$E_x$', transform=ax1.transAxes, fontsize=18)
ax2.text(0.85, 0.9, r'$E_y$', transform=ax2.transAxes, fontsize=18)
ax3.text(0.85, 0.9, r'$E_z$', transform=ax3.transAxes, fontsize=18)
ax4.text(0.85, 0.9, r'$B_x$', transform=ax4.transAxes, fontsize=18)
ax5.text(0.85, 0.9, r'$B_y$', transform=ax5.transAxes, fontsize=18)
ax6.text(0.85, 0.9, r'$B_z$', transform=ax6.transAxes, fontsize=18)
cb1 = fig.colorbar(p3, ax=ax3)
cb1.set_label(r'[V/m]',rotation=0, position=[0.0,1.0])
cb2 = fig.colorbar(p6, ax=ax6)
cb2.set_label(r'[T]',rotation=0, position=[0.0,1.0])
        
plt.show()

