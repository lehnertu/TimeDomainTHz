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

#include <Eigen/Dense>

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

    // **************************************
    // load the input field from a file
    // **************************************

    std::cout << std::endl << "=== Source Screen ===" << std::endl;
    Screen *source = new Screen(infile);
    source->init();
    // print report
    source->writeReport(&cout);

    Vector avg_normal = Vector(0.0,0.0,0.0);
    Vector avg_xi = Vector(0.0,0.0,0.0);
    Vector avg_eta = Vector(0.0,0.0,0.0);
    for (int ip=0; ip<source->get_Np(); ip++)
    {
        avg_normal += source->normal[ip];
        avg_xi += source->xi[ip];
        avg_eta += source->eta[ip];
    };
    std::cout << "n =   (" << avg_normal.x << ", " << avg_normal.y << ", " << avg_normal.z << ")" << std::endl;
    std::cout << "xi =  (" << avg_xi.x << ", " << avg_xi.y << ", " << avg_xi.z << ")" << std::endl;
    std::cout << "eta = (" << avg_eta.x << ", " << avg_eta.y << ", " << avg_eta.z << ")" << std::endl;

    // **************************************
    // test the computation of derivatives
    // **************************************
    // temporarily we define most variables of the screen as public

    // test the EIGEN library
    Eigen::MatrixXf A = Eigen::MatrixXf::Random(3, 2);
    cout << "Here is the matrix A:\n" << A << endl;
    Eigen::VectorXf b = Eigen::VectorXf::Random(3);
    cout << "Here is the right hand side b:\n" << b << endl;
    cout << "The least-squares solution is:\n"
        << A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b) << endl;

    std::cout << "computing the derivatives of the fields ..." << std::endl;
    // record the start time
    timespec start_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    
    int index = 1190;
    int *nbh = new int[20];
    int count = source->get_Neighbourhood(index,nbh);
    std::cout << "neighbourhood size = " << count << std::endl;
    std::cout << "[";
    for (int i=0; i<count; i++) std::cout << nbh[i] << " ";
    std::cout << "]" << std::endl;
    
    // an array of 20 pointers to field traces in the local coordinate system
    FieldTrace* local_trace[20];
    for (int i=0; i<count; i++)
    {
        local_trace[i] = new FieldTrace(source->A[nbh[i]]);
        local_trace[i]->transform(source->xi[index],source->eta[index],source->normal[index]);
    };
    
    for (int i=0; i<count; i++) delete local_trace[i];
    delete nbh;

    // record the finish time
    timespec stop_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    double elapsed = stop_time.tv_sec-start_time.tv_sec +
        1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    std::cout << "time elapsed during computation : " << elapsed << " s" << std::endl;
    
    // **************************************
    // load the target geometry from a file
    // **************************************

    std::cout << std::endl << "=== Target Screen ===" << std::endl;
    Screen *target = new Screen(outfile);
    target->init();
    // print report
    target->writeReport(&cout);
    
    avg_normal = Vector(0.0,0.0,0.0);
    avg_xi = Vector(0.0,0.0,0.0);
    avg_eta = Vector(0.0,0.0,0.0);
    for (int ip=0; ip<target->get_Np(); ip++)
    {
        avg_normal += target->normal[ip];
        avg_xi += target->xi[ip];
        avg_eta += target->eta[ip];
    };
    std::cout << "n =   (" << avg_normal.x << ", " << avg_normal.y << ", " << avg_normal.z << ")" << std::endl;
    std::cout << "xi =  (" << avg_xi.x << ", " << avg_xi.y << ", " << avg_xi.z << ")" << std::endl;
    std::cout << "eta = (" << avg_eta.x << ", " << avg_eta.y << ", " << avg_eta.z << ")" << std::endl;

    double min_t0 = 1.0e30;
    double max_t0 = -1.0e30;
    for (int i=0; i<target->get_Np(); i++)
    {
        double t0 = target->get_t0(i);
        if (t0>max_t0) max_t0=t0;
        if (t0<min_t0) min_t0=t0;
    }
    std::cout << "target start timing range (" << min_t0 << ", " << max_t0 << ")" << std::endl;
    
    delete source;
    delete target;
    
    return 0;
}
