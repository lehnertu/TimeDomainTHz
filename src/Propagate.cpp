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
    \brief Time-Domain THz general (screen-to-screen) propagation tool.
    
    @author Ulf Lehnert
    @date 21.11.2019
    @file Propagate.cpp
    
    The input beam is read from a HDF5 file defining the field traces on a meshed screen.
    The geometry of the output screen is read from a second HDF5 file.
    From this second file only the geometry information is used, it will be
    rewritten with the computed output fields.

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
    std::cout << "Propagate Screen to Screen" << std::endl << std::endl;
    
    if (argc<3)
    {
        std::cout << "Usage: Propagate infile outfile" << std::endl;
        return 1;
    }
    std::string infile(argv[1]);
    std::string outfile(argv[2]);

    // load the input field from a file
    Screen *source = new Screen(infile);
    // print report
    source->writeReport(&cout);
    source->writeTraceReport(&cout, 0);

    delete source;
    
    return 0;
}
