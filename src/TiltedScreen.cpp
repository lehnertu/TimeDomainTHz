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
    \brief Propagate Input beam to a screen which is tilted to the direction of incidence

    @author Ulf Lehnert
    @date 21.11.2019
    @file Gaussian-TiltedScreen.cpp
    
    The time-domain source field is loaded from an HDF5 file.
    The beam is propagating in positive z-direction.
    Polarization of the electric field is horizontal (x-direction).
    
    The orientation vectors of the screen are positive x (right)
    and negative y direction (up) so the normal vector is negative z-direction.
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
    std::cout << "Propagate to Tilted Screen" << std::endl << std::endl;

    // load the input field from a file
    if (argc<2)
    {
        std::cout << "Error - file name parameter required" << std::endl;
        return 1;
    }
    std::string filename(argv[1]);
    Screen *source = new Screen(filename);
    // print report
    source->writeReport(&cout);
    // compute the derivatives of the source fields
    std::cout << "computing derivatives of source field..." << std::endl;
    source->computeDerivatives();
    std::cout << std::endl;
    // write the screen data to file
    // with modified sources we write the derivatives insted
    // Screen::bufferArray() is changed to deliver dx_A instead of A
    // source->writeFieldHDF5("Gaussian_51_dx.h5");

    // define the target screen
    double distance = 1.25751;
    int Nx = 51;
    int Ny = 51;
    FieldTrace source_trace = source->get_Trace(0,0);
    int Nt = source_trace.get_N()+400;
    Screen *target = new Screen(
        Nx, Ny, Nt,
        Vector(0.002,0.0,0.0),
        Vector(0.0,-0.002,0.0),
        source->get_Center() + Vector(0.0,0.0,distance) );
    target->writeReport(&cout);

    // compute the source beam propagated to the target screen
    std::cout << "computing propagated source fields ..." << std::endl;
    double target_t0 = source->get_t0()+distance/SpeedOfLight;
    for (int ix=0; ix<Nx; ix++)
    {
        for (int iy=0; iy<Ny; iy++)
        {
            Vector pos = target->get_point(ix,iy);
            FieldTrace target_trace = source->propagation(pos, target_t0, Nt);
            target->set_Trace(ix, iy, target_trace);
            std::cout << iy << " ";
        };
        std::cout << ix << std::endl;
    };
    std::cout << "done." << std::endl;
    std::cout << std::endl;

    // write the target screen data to file
    target->writeFieldHDF5("Mirror_25.h5");
    
    delete source;
    delete target;
    
    return 0;
}
