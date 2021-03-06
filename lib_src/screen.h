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

/*  Grid vectors are not orthogonal to each other -> throw an exception */
class Screen_NonOrthogonal { };
/*  Memory allocation may go wrong -> throw an exception */
class Screen_MemoryAllocationError { };
/*  Assignments of field traces may go wrong -> throw an exception */
class Screen_IndexOutOfRange { };
/*  Writing an HDF5 file may go wrong -> throw an exception */
class Screen_FileWriteError { };
/*  Reading an HDF5 file may go wrong -> throw an exception */
class Screen_FileReadError { };
/*  further propagate a FieldTrace_Zero exception */
struct Screen_Zero_exception {
    int ix;
    int iy;
    Screen_Zero_exception(int x, int y) {ix=x; iy=y;};
};

/*!
 * \class Screen
 * \brief rectangular array of electromagnetic field traces
 * @author Ulf Lehnert
 * @date 18.11.2019
 * 
 * This class holds electromagnetic field traces on a rectangular grid.
 * All traces are supposed to have the same temporal definitions.
 */
class Screen
{

public:

    /*! Default constructor:<br>
     *  @param Nx number of grid cells in xVec direction
     *  @param Nx number of grid cells in xVec direction
     *  @param Nt number of samples in time
     *  @param xVec defines x-axis of grid spacing
     *  @param yVec defines y-axis of grid spacing
     *  @param Center center point of the screen
     *
     *  xVec and yVec are restricted to be orthogonal.
     */
    Screen(
        int Nx, int Ny, int Nt,
        Vector xVec, Vector yVec,
        Vector center
        );
        
    /*! Construct a screen object from an HDF5 file */
    Screen(std::string filename);
    
    /*! Default destructor:<br>
     */
    ~Screen();

    /*! Getter routines for members
     *  no setter routines are provided,
     *  these values are not expected to change once constructed */
    int get_Nx() { return Nx; }
    int get_Ny() { return Ny; }
    int get_Nt() { return get_trace(0,0).get_N(); }
    double get_dt() { return get_trace(0,0).get_dt(); }
    Vector get_Center()  { return Center; }
    Vector get_Normal()  { return Normal; }
    double get_dA()  { return dA; }
    Vector get_xVec()  { return xVec; }
    Vector get_yVec()  { return yVec; }
    
    /*! Get the position of one grid cell */
    Vector get_point(int ix, int iy); 
    
    /*! Get the field trace data for one of the grid cells */
    FieldTrace get_trace(int ix, int iy);

    /*! Set the field trace data for one of the grid cells */
    void set_trace(int ix, int iy, FieldTrace trace);
    
    /*! Total energy of radiation falling on the screen */
    double totalEnergy();
    
    /*! Write a report of the screen geometry and parameters
     *  including some summary field data
     *  onto an output stream
     */
    void writeReport(std::ostream *st);
    
    /*! Write a report of a single field trace
     *  onto an output stream
     */
    void writeTraceReport(std::ostream *st, int ix, int iy);

    /*! Determine the size of a buffer needed to hold the data array
     *  @return number of doubles
     */
    int getBufferSize() { return Nx*Ny*A[0][0].get_N()*6; }
    
    /*! Copy the data array into a given buffer.
     *  The size of the buffer has to be as determined by getBufferSize().
     */
    void bufferArray(double *buffer);

    /*! Compute the spatial and temporal derivatives of the fields.
     *  internally parallelized using Open-MP
     *  @TODO The computation is flawed as it assumes orthogonality of xVec and yVec
     *  when computing the normal derivatives.
     *  Signs of the drivatives have been checked to match those of the python code.
     *  This may be wrong, in particular for dy_A (because yVec is -ey in the test case)
     *  and dn_A should be dz_A but is set to -dz_A.
     *  @TODO For the propagation only the temporal and normal derivatives are needed.
     *  It is not necessary to store the transverse derivatives.
     */    
    void computeDerivatives();

    /*! Compute the fields propagated to a given point in space.
     *  It is assumed that the normal vector of the source area points outward from the volume of interest.
     *  The time sampling of the target trace is retained as given.
     *  The fields are properly retarded at the target point according to its distance from the source.
     *  @param target_pos observation position for the target trace
     *  @param target_trace field trace to be overwritten with the result
     */
    void propagate_to(Vector target_pos, FieldTrace *target_trace);

    /*! Write the screen data to an HDF5 file */
    void writeFieldHDF5(std::string filename);

private:

    /*! Compute the field derivative along the xVec direction.
     *  This is used when computing the normal derivatives.
     */
    FieldTrace dx_A(int ix, int iy);

    /*! Compute the field derivative along the yVec direction.
     *  This is used when computing the normal derivatives.
     */
    FieldTrace dy_A(int ix, int iy);

private:

    /*! number of grid cells in each direction */
    int Nx, Ny;
    
    /*! grid vectors in each direction */
    Vector xVec, yVec;
    
    /*! area element of a single grid point */
    double dA;
    
    /*! normal vector of the screen pointing to the incidence (backward) side */
    Vector Normal;
    
    /*! center of the screen */
    Vector Center;
    
    /*! array [Nx][Ny] of FieldTrace
     *  outer (first) index is Nx
    */
    std::vector< std::vector<FieldTrace> > A;

    /*! partial derivative of the fields along the normal direction
     *  data are only valid after calling computeDerivatives()
     */
    std::vector< std::vector<FieldTrace> > dn_A;

    /*! time derivative of the fields
     *  data are only valid after calling computeDerivatives()
     */
    std::vector< std::vector<FieldTrace> > dt_A;
};

