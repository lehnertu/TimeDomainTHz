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
#include <time.h>
#include <omp.h>

#include "global.h"
#include "fields.h"
#include "screen.h"

int main(int argc, char* argv[])
{
    std::cout << std::endl;
    std::cout << "Propagate to Screen" << std::endl << std::endl;
    
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
    double distance = 2.0 * 1.25751;
    int Nx = 51;
    int Ny = 51;
    FieldTrace source_trace = source->get_Trace(0,0);
    int Nt = source_trace.get_N()+4000;
    Screen *target = new Screen(
        Nx, Ny, Nt,
        Vector(0.003,0.0,0.0),
        Vector(0.0,-0.003,0.0),
        source->get_Center() + Vector(0.0,0.0,distance) );
    target->writeReport(&cout);

    // compute the source beam propagated to the target screen
    double target_t0 = source->get_tCenter()+distance/SpeedOfLight-source->get_dt()*(double)Nt/2.0;
    time_t start_t;
    time(&start_t);
    time_t print_time = start_t;
    time_t now;
    // count the target screen points already computed
    int counter = 0;
    // parallel domain
    // optional parameter :  num_threads(32)
    // all variables declared outside this block are shared, e.g. source, target
    // all variables declared inside this block are private, e.g. pos, target_trace
    #pragma omp parallel shared(counter)
    {
        #pragma omp single
        {
            std::cout << "computing propagated source fields on " << omp_get_num_threads() << " parallel threads" << std::endl;
        }
        #pragma omp for
        for (int i=0; i<Nx*Ny; i++)
        {
            // we comprise the two loops over ix and iy into one parallel loop
            int ix = i/Ny;
            int iy = i - ix*Ny;
            /* pragma omp critical
            {
                std::cout << "ix = " << ix << "  iy = " << iy << std::endl;
            };
            */
            Vector pos = target->get_point(ix,iy);
            FieldTrace target_trace = source->propagation(pos, target_t0, Nt);
            target->set_Trace(ix, iy, target_trace);
            // all threads increment the counter
            #pragma omp atomic
            counter++;
            // only the master thread keeps the time
            if (omp_get_thread_num()==0)
            {
                time(&now);
                if (difftime(now, print_time)>10.0)
                {
                    print_time = now;
                    std::cout << "completed " << 100.0*(double)counter/(double)(Nx*Ny) << "%";
                    std::cout << "   after " << difftime(now, start_t) << " seconds" << std::endl;
                };
            };
        };
    };
    // end parallel domain
    time(&now);
    std::cout << "completed after " << difftime(now, start_t) << "seconds." << std::endl;
    std::cout << std::endl;
    target->writeReport(&cout);

    // write the target screen data to file
    target->writeFieldHDF5("Gaussian_25_Straight_51_C.h5");
    
    delete source;
    delete target;
    
    return 0;
}