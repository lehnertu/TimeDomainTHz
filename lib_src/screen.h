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

#pragma once

#include <iostream>

#include "vector.h"
#include "fields.h"

/*  Memory allocation may go wrong -> throw an exception */
class Screen_MemoryAllocationError { };
/*  Assignments of field traces may go wrong -> throw an exception */
class Screen_IndexOutOfRange { };
/*  Writing an HDF5 file may go wrong -> throw an exception */
class Screen_FileWriteError { };
/*  Reading an HDF5 file may go wrong -> throw an exception */
class Screen_FileReadError { };

/*!
 * \class Screen
 * \brief rectangular array of electromagnetic field traces
 * @author Ulf Lehnert
 * @date 18.11.2019
 * 
 * This class holds electromagnetic field traces on an unstructured mesh.
 * The mesh consists of triangles, the center points are the positions at which
 * the fields are defined. All traces are have a common length and time-step but
 * can have differing starting times.
 */
class Screen
{

public:

    /*! Default constructor:<br>
     *  @param Ncp number of the corner points of the mesh.
     *  @param Np number of the mesh cells and field traces
     *  @param Nt number of samples in time
     */
    Screen(
        int Ncp, int Np, int Nt,
        double dt
        );

    /*! Construct a screen object from an HDF5 file */
    Screen(std::string filename);
    
    /*! Default destructor:<br>
     */
    ~Screen();

    /*! Getter routines for members
     *  no setter routines are provided,
     *  these values are not expected to change once constructed */
    int get_Ncp() { return Ncp; }
    int get_Np() { return Np; }
    int get_Nt() { return Nt; }
    double get_dt() { return dt; }

    /*! Get the position of one grid cell */
    Vector get_point(int ip) { return field_points[ip]; }

    /*! Write a report of the screen geometry and parameters
     *  including some summary field data
     *  onto an output stream
     */
    void writeReport(std::ostream *st);
    
    /*! Write a report of a single field trace
     *  onto an output stream
     */
    void writeTraceReport(std::ostream *st, int ip);
        
private:

    /*! number of corner points of the grid cells */
    int Ncp;
    
    /*! number of the grid cells / field observation positions */
    int Np;

    /*! number of the time steps per trace */
    int Nt;
    
    /*! time-step value */
    double dt;

    /*! starting time of all the traces */
    std::vector<double> t0;

    /*! the coordinates of the mesh corner points */
    std::vector<Vector> triangle_points;
    
    /*! the coordinates of the mesh center points where fields are defined */
    std::vector<Vector> field_points;

    /*! field traces for every mesh cell */
    std::vector<FieldTrace> A;

};

