#!/usr/bin/env python3
# coding=UTF-8

import sys, time
import os
import psutil
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Circle

pid = os.getpid()
print("PID = ",pid)
libmem = psutil.Process(pid).memory_info().rss/2.0**20
print('libraries memory use: %.2f MB' % libmem)

# magnetic field constant in N/A²
mu0 = 4*np.pi*1e-7

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the screen output HDF5 file')
parser.add_argument('-xy', help="indeces of plot point", dest="xy", type=int, nargs=2)

print()
args = parser.parse_args()

radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print("file not found")
  sys.exit()

# Open the file for reading
print("reading %s" % radfile)
hdf = h5py.File(radfile, "r")
print(hdf)
# Get the groups
pos = hdf['ObservationPosition']
Nx = pos.attrs.get('Nx')
Ny = pos.attrs.get('Ny')
xcenter = (Nx-1)//2
ycenter = (Ny-1)//2
print("Nx=%d Ny=%d" % (Nx,Ny))
print(pos)
time = hdf['ObservationTime']
nots = time.attrs.get('Nt')
dt = time.attrs.get('dt')
print("Nt=%d dt=%g" % (nots,dt))
print(time)
field = hdf['ElMagField']
print(field)
pos = np.array(pos)
t0 = np.array(time)
# a = np.array(field)
# direct slicing to reduce the memory requirements
# this should never store the complete array
onaxis = np.array(field[xcenter,ycenter,:,:])
if args.xy != None:
    xi = args.xy[0]
    yi = args.xy[1]
    offaxis = np.array(field[xi,yi,:,:])
hdf.close()
print("file closed.")

print('data memory use: %.2f MB' % (psutil.Process(pid).memory_info().rss/2.0**20-libmem))
print()

print("center = (%d, %d)" % (xcenter,ycenter))
centerposition = pos[xcenter][ycenter]
print("center position = %s" % centerposition)
print("start time = %.4f ns" %(1e9*t0[xcenter][ycenter]))

Ex = onaxis[:,0]
Ey = onaxis[:,1]
Ez = onaxis[:,2]
Bx = onaxis[:,3]
By = onaxis[:,4]
Bz = onaxis[:,5]

EVec = np.array([Ex, Ey, Ez]).transpose()
BVec = np.array([Bx, By, Bz]).transpose()
# Poynting vector in V/m * (N/(A m)) / (N/A²) = W/m²
SVec = np.cross(EVec, BVec) / mu0
t = 1e9*np.linspace(t0[xcenter][ycenter],t0[xcenter][ycenter]+(nots-1)*dt,nots)
print("on axis energy flow density = %s µJ/m²" % (1e6*SVec.sum(axis=0)*dt))

if args.xy != None:

    print("index = (%d, %d)" % (xi,yi))
    position = pos[xi][yi]
    print("off-axis position = %s" % position)
    print("start time = %.4f ns" %(1e9*t0[xi][yi]))

    Ex = offaxis[:,0]
    Ey = offaxis[:,1]
    Ez = offaxis[:,2]
    Bx = offaxis[:,3]
    By = offaxis[:,4]
    Bz = offaxis[:,5]

    EVec = np.array([Ex, Ey, Ez]).transpose()
    BVec = np.array([Bx, By, Bz]).transpose()
    # Poynting vector in V/m * (N/(A m)) / (N/A²) = W/m²
    SVec = np.cross(EVec, BVec) / mu0
    t = 1e9*np.linspace(t0[xi][yi],t0[xi][yi]+(nots-1)*dt,nots)
    print("off axis energy flow density = %s µJ/m²" % (1e6*SVec.sum(axis=0)*dt))

# prepage the plot

left, width = 0.15, 0.80
rect1 = [left, 0.55, width, 0.40]  #left, bottom, width, height
rect2 = [left, 0.08, width, 0.40]
fig = plt.figure(1,figsize=(12,9))

ax1 = fig.add_axes(rect1)
ax4 = fig.add_axes(rect2, sharex=ax1)

# plot the time-trace of the fields

l1 = ax1.plot(t, Ex, "r-", label=r'$E_x$')
l2 = ax1.plot(t, Ey, "b-", label=r'$E_y$')
l3 = ax1.plot(t, Ez, "g-", label=r'$E_z$')

ax1.set_ylabel(r'$E$ [V/m]')
lines = l1 + l2 + l3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc='upper right')
for label in ax1.get_xticklabels():
    label.set_visible(False)
ax1.grid(True)

l4 = ax4.plot(t, Bx, "r-", label=r'$B_x$')
l5 = ax4.plot(t, By, "b-", label=r'$B_y$')
l6 = ax4.plot(t, Bz, "g-", label=r'$B_z$')

ax4.set_ylabel(r'$B$ [T]')
ax4.set_xlabel(r't [ns]')
lines = l4 + l5 +l6
labels = [l.get_label() for l in lines]
ax4.legend(lines,labels,loc='upper right')
ax4.grid(True)

plt.show()
