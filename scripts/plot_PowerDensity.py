#!/usr/bin/env python3
# coding=UTF-8

import sys, time
import os.path
import argparse
import numpy as np
from TimeDomainField import *
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Circle

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the input file')

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
print("dt=%g Nt=%d" % (Source.dt, Source.Nt))
print()

xc, yc = Source.center_index()
center = Source.center()
print("center = (%d, %d)" % (xc,yc))
print("position = %s" % center)
print()

print("dxVec = %s m   dX = %s m" % (Source.dxVec, Source.dX))
print("dyVec = %s m   dY = %s m" % (Source.dyVec, Source.dY))
print("n = %s" % Source.normal)
print("pixel size = (%9.6f, %9.6f) %6.3g m²" % (Source.dX,Source.dY,Source.dX*Source.dY) )
print()

print('energy flow (Poynting) vector on axis = %s J/m²' % Source.poyntingVector(xc, yc))
print('normal energy flow density at center = %s J/m²' % Source.powerDensity(xc, yc))
print('total energy = %s  µJ' % (1e6*Source.totalPower()))

# plot figure with power density on screen
# density array has to be transposed
Pn = np.zeros([Source.Ny, Source.Nx])
x = np.zeros([Source.Nx])
for ix in range(Source.Nx):
    x[ix] = np.dot(Source.pos[ix,yc]-center,Source.dxVec/Source.dX)
y = np.zeros(Source.Ny)
for iy in range(Source.Ny):
    y[iy] = np.dot(Source.pos[xc,iy]-center,Source.dyVec/Source.dY)
for ix in range(Source.Nx):
    for iy in range(Source.Ny):
        Pn[iy,ix] = Source.powerDensity(ix,iy)

fig1 = plt.figure(1,figsize=(12,9))
ax1 = fig1.add_subplot(111)
plt.contourf(x, y, Pn, 15, cmap='CMRmap')
plt.title('energy flow density [J/m$^2$]')
plt.xlabel('x /m')
plt.ylabel('y /m')
cb=plt.colorbar()
cb.set_label(r'energy density [J/m$^2$]')

plt.show()

