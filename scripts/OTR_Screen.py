#!/usr/bin/env python
# coding=UTF-8

import sys, time
import os.path
import numpy as np
from scipy import constants
import h5py
import matplotlib.pyplot as plt

"""
We compute the transient fields of a single point charge
passing through the origin of a screen at time zero.

The particle shall have relativistic velocity (gamma given)
and move in positive z direction. The screen normal points
in negative z direction with the transverse axes being
x-left and y-up.
"""

gamma = 3.0
betagamma = np.sqrt(gamma*gamma-1.0)
beta = np.sqrt(1.0-1.0/(gamma*gamma))

# Lorentz transformation 4-tensor
Lorentz_Z = np.array(
    [[gamma,0.0,0.0,-betagamma],
     [0.0,1.0,0.0,0.0],
     [0.0,0.0,1.0,0.0],
     [-betagamma,0.0,0.0,gamma]])
Inv_Lorentz_Z = np.linalg.inv(Lorentz_Z)
print(Lorentz_Z)
print(Inv_Lorentz_Z)

# particle position in lab frame at a given time
# (c*t, x, y, z)
def Lab_X (t):
    return np.array([constants.c*t,0.0,0.0,beta*constants.c*t]);

# particle position transformed to co-moving frame
# (c*t', x', y', z')
def Mov_X (X):
    return Lorentz_Z.dot(X)

# (static) electric field of the particle in the co-moving frame
# given as 4-tensor
# parameter is the position in the co-moving frame
def Mov_F (X):
    f = -constants.e / (4.0*np.pi*constants.epsilon_0)
    x = X[1]
    y = X[2]
    z = X[3]
    R2 = x*x + y*y + z*z
    Ex_c = f*x/R2/constants.c
    Ey_c = f*y/R2/constants.c
    Ez_c = f*z/R2/constants.c
    return np.array(
        [[0.0,-Ex_c,-Ey_c,Ez_c],
        [Ex_c,0.0,0.0,0.0],
        [Ey_c,0.0,0.0,0.0],
        [Ez_c,0.0,0.0,0.0]])

# grid parameters
Nz = 100
Nx = 100
span_Z = 0.001
span_X = 0.001
# position on screen
z = np.linspace(-span_Z/2.0,span_Z/2.0,num=Nz)
x = np.linspace(-span_X/2.0,span_X/2.0,num=Nx)

A = np.zeros((Nz,Nx))

for iz in range(Nz):
    for ix in range(Nx):
        obs_pos = np.array([0.0,x[ix],0.0,z[iz]])
        local_pos = Mov_X(obs_pos)
        mov_F = Mov_F(local_pos)
        obs_F = Inv_Lorentz_Z.dot(Inv_Lorentz_Z.dot(mov_F))
        Ex = obs_F[1,0]*constants.c
        Ey = obs_F[2,0]*constants.c
        Ez = obs_F[3,0]*constants.c
        # A[iz,ix] = np.sqrt(Ex*Ex + Ey*Ey + Ez*Ez)
        A[iz,ix] = np.sqrt(Ex*Ex)
        # A[iz,ix] = np.sqrt(Ez*Ez)

# plot data
fig1 = plt.figure(1,figsize=(12,9))
ax1 = fig1.add_subplot(111)
plt.contourf(z, x, np.transpose(A), 15, cmap='CMRmap')
plt.title('electric field magnitude [V/m]')
plt.xlabel('z /m')
plt.ylabel('x /m')
cb=plt.colorbar()
cb.set_label(r'field [V/m]')

plt.show()

"""
if __name__ == '__main__':

    # write HDF5 file

    hf = h5py.File('Forward_OTR_25.h5', 'w')
    h5p = hf.create_dataset('ObservationPosition', data=pos)
    h5p.attrs['Nx'] = Nx
    h5p.attrs['Ny'] = Ny
    h5f = hf.create_dataset('ElMagField', data=A)
    h5f.attrs['t0'] = 0.0
    h5f.attrs['dt'] = dt
    h5f.attrs['NOTS'] = NOTS
    hf.close()

"""
