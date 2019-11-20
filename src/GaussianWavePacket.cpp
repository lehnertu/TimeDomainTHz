/*=========================================================================
 * 
 *  Program:   TimeDomainTHz - THz radiation transport in time-domain
 *
 *  Copyright (c) 2019 U. Lehnert
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * =========================================================================*/

/*!
    \brief Gaussian beam with gaussian time dependence

    @author Ulf Lehnert
    @date 18.11.2019
    @file GaussianWavePacket.cpp
    
    Create a time-domain field of a Gaussian wave-packet on a screen.
    Save as an HDF5 file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "global.h"
#include "fields.h"
#include "screen.h"

int main(int argc, char* argv[])
{
    std::cout << std::endl;
    std::cout << "Gaussian Wave Packet" << std::endl << std::endl;

    // define the wave parameters
    double tau = 4.0e-12;
    double f = 0.3e12;
    double w0 = 0.02;
    double lam = SpeedOfLight/f;
    std::cout <<  "lambda = " << lam*1e6 << " µm" << std::endl;
    double zR = Pi*w0*w0/lam;
    std::cout <<  "w0 = " << 1e3*w0 << " mm" << std::endl;
    std::cout <<  "Rayleigh range = " << zR << " m" << std::endl << std::endl;
    
    // define the peak field
    // electric field amplitude in V/m
    // magnetic field amplitude in T
    ElMagField peak = ElMagField(Vector(1.0,0.0,0.0),Vector(0.0,1.0/SpeedOfLight,0.0)) * 2.0e7;

    // compute the fields on axis including all time dependencies
    // this will be the same for all grid points, only field strength varies
    FieldTrace *on_axis = new FieldTrace(0.0, 1.0e-13, 400);
    for (int i=0; i<on_axis->get_N(); i++)
    {
        double t = on_axis->get_time(i);
        double osc = cos(2.0*Pi*f*t) * exp(-(t*t)/(2.0*tau*tau));
        on_axis->set(i, peak*osc);
    };
    Vector dPdA = on_axis->Poynting();
    std::cout << "Power density on axis = ";
    std::cout << "(" << dPdA.x << ", " << dPdA.y << ", " << dPdA.z << ")";
    std::cout << " J/m²" << std::endl;
    
    // setup the geometry of the screen
    int Nx = 51;
    int Ny = 51;
    Screen *scr = new Screen( Nx, Ny, 400,
        Vector(0.05,0.0,0.0), Vector(0.0,0.05,0.0),
        Vector(0.0,0.0,0.0) );
    
    // set the field traces for all grid points of the screen
    Vector center = scr->get_Center();
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            Vector pos = scr->get_point(ix,iy);
            double r = (pos-center).norm();
            double radint = exp(-r*r/(w0*w0));
            scr->set_Trace(ix,iy, (*on_axis)*radint);
        };
    std::cout << "Energy incident on screen = " << scr->totalEnergy() << " J" << std::endl;
    
    // write the screen data to file
    /*
    hf = h5py.File('Gaussian_25.h5', 'w')
    h5p = hf.create_dataset('ObservationPosition', data=pos)
    h5p.attrs['Nx'] = Nx
    h5p.attrs['Ny'] = Ny
    h5f = hf.create_dataset('ElMagField', data=A)
    h5f.attrs['t0'] = 0.0
    h5f.attrs['dt'] = dt
    h5f.attrs['NOTS'] = NOTS
    hf.close()
    */

    return 0;
}
