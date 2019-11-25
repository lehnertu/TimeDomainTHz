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

    // define the target screen
    double distance = 1.25751;
    int Nx = 25;
    int Ny = 25;
    FieldTrace source_trace = source->get_Trace(0,0);
    int Nt = source_trace.get_N()+400;
    Screen *target = new Screen(
        Nx, Ny, Nt,
        Vector(0.005,0.0,0.0),
        Vector(0.0,-0.005,0.0),
        source->get_Center() + Vector(0.0,0.0,distance) );
    target->writeReport(&cout);

    delete source;
    return 0;
}
