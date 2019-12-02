import numpy as np
from scipy import constants
import h5py

# *******************************************
# class definition for time-domain fields
# *******************************************

class TimeDomainField:

    def __init__(self, nx=0, ny=0, nt=0):
        self.Nx = nx
        self.Ny = ny
        self.Nt = nt
        if (self.Nx>0) and (self.Ny>0) and (self.Nt>0):
            self.A = np.zeros((self.Nx, self.Ny, self.Nt, 6))
        
    # read all information from an HDF5 file
    def read(self, filename):
        hdf = h5py.File(filename, "r")
        p = hdf['ObservationPosition']
        self.Nx = p.attrs.get('Nx')
        self.Ny = p.attrs.get('Ny')
        self.pos = np.array(p)
        self.dxVec = self.pos[1,0]-self.pos[0,0]
        self.dX = np.linalg.norm(self.dxVec)
        self.dyVec = self.pos[0,1]-self.pos[0,0]
        self.dY = np.linalg.norm(self.dyVec)
        # this is the backward normal assuming z-direction beam propagation
        # x-left  y-up  z-forward
        self.normal = np.cross(self.dxVec,self.dyVec)
        self.normal = self.normal / np.linalg.norm(self.normal)
        t = hdf['ObservationTime']
        self.dt = t.attrs.get('dt')
        self.Nt = t.attrs.get('Nt')
        self.t0 = np.array(t)
        f = hdf['ElMagField']
        self.A = np.array(f)
        hdf.close()

    # write all information to an HDF5 file
    def write(self, filename):
        hf = h5py.File(filename, 'w')
        h5p = hf.create_dataset('ObservationPosition', data=self.pos)
        h5p.attrs['Nx'] = self.Nx
        h5p.attrs['Ny'] = self.Ny
        h5p = hf.create_dataset('ObservationTime', data=self.t0)
        h5p.attrs['Nt'] = self.Nt
        h5p.attrs['dt'] = self.dt
        h5f = hf.create_dataset('ElMagField', data=self.A)
        hf.close()

    # return the center of the field distribution
    def center_index(self):
        ixc = (self.Nx-1)//2
        iyc = (self.Ny-1)//2
        return ixc, iyc
    def center(self):
        ixc, iyc = self.center_index()
        return self.pos[ixc,iyc]
    
    # compute the (vectorial) power flow density at one point on the grid
    def poyntingVector(self, ix, iy):
        trace = self.A[ix,iy]
        EVec = trace[:,0:3]
        BVec = trace[:,3:6]
        SVec = np.cross(EVec, BVec) / constants.mu_0
        P = (SVec.sum(axis=0))*self.dt
        return P
        
    # compute the power density (per grid area) at one point on the grid
    # component in the opposite direction of the normal vector
    def powerDensity(self, ix, iy):
        return -np.dot(self.poyntingVector(ix,iy),self.normal)

    # compute the total power of the field
    def totalPower(self):
        total = 0.0
        for ix in range(self.Nx):
            for iy in range(self.Ny):
                total += self.powerDensity(ix,iy)
        return total*self.dX*self.dY

