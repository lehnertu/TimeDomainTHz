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

// TODO: all index out-of-range possibilities have to be checked

/*  Writing an HDF5 file may go wrong -> throw an exception */
class Screen_FileWriteError { };
/*  Reading an HDF5 file may go wrong -> throw an exception */
class Screen_FileReadError { };
/*  If normals are not sifficiently aligned -> throw an exception */
class Screen_NormalsAlignmentError { };
/*  Memory allocation may go wrong -> throw an exception */
class Screen_MemoryAllocationError { };

/*! data type for tringle references */
struct tri_ref {
    int p1;
    int p2;
    int p3;
};

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
    
    /*! This method computes a number of internal data
     *  like cell areas, normals and performs some sanity checks.
     *  It should be called whenever a new screen object has been constructed.
     *
     *  A local coordinate system is set up for every mesh cell.
     *  The direction vectors (xi, eta, n) form a right-handed cartesic system.
     */
    void init();
    
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

    /*! Get the start time of one trace */
    double get_t0(int ip) { return t0[ip]; }

    /*! Get the cell-local coordinate system */
    Vector get_xi(int ip) { return xi[ip]; }
    Vector get_eta(int ip) { return eta[ip]; }
    Vector get_normal(int ip) { return normal[ip]; }

    /*! Get the area of one cell */
    double get_area(int ip) { return area[ip]; }

    /*! Get the fields of one cell */
    FieldTrace get_trace(int ip) { return *A[ip]; }

    /*! Set the fields of one cell */
    void set_trace(int ip, FieldTrace trace) { *A[ip]=trace; };
    
    /*! Determine the neighbourhood of a cell with given index.
     *  These are all cells that share at least one common point with the indexed cell.
     *  A buffer of sufficient size to hold the cell indices (up to 20) has to be provided.
     *  Return value is the number of cells in the neighbourhood (including the indexed cell).
     */
    int get_Neighbourhood(int index, int *buffer);
    
    /*! Get the energy flow density vector of one trace
     *  integrated over time
     */
    Vector Poynting(int ip) { return A[ip]->Poynting(); };

    /*! Compute the total electromagnetic energy flowing through the screen
     *  Energy flow vectors opposite the normal vector are counted positive.
     */
    double totalEnergy();
    
    /*! Write a report of the screen geometry and parameters
     *  including some summary field data
     *  onto an output stream
     */
    void writeReport(std::ostream *st);
    
    /*! Write a report of a single field trace
     *  onto an output stream
     */
    void writeTraceReport(std::ostream *st, int ip);

    /*! Write the screen data to an HDF5 file */
    void writeFieldHDF5(std::string filename);
    
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
    
    /*! lists of references which points belong to the triangles */
    std::vector<tri_ref> triangles;
    
    /*! the coordinates of the mesh center points where fields are defined */
    std::vector<Vector> field_points;

    /*! the local coordinate system of the cells - computed by init() */
    std::vector<Vector> xi;
    std::vector<Vector> eta;
    std::vector<Vector> normal;

    /*! the cell areas - computed by init() */
    std::vector<double> area;
    
    /*! the total screen area - computed by init() */
    double total_area;

    /*! pointers to the field traces for every mesh cell */
    std::vector<FieldTrace*> A;

};

