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

#include "vector.h"
#include "fields.h"

/*  Assignments of field traces may go wrong -> throw an exception */
class Screen_IndexOutOfRange { };

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
     *  @param yVec defines y-axis of grid spacing (needs not to be perpendicular to xVec)
     *  @param Center center point of the screen
     */
    Screen(
        int Nx, int Ny, int Nt,
        Vector xVec, Vector yVec,
        Vector center
        );

    /*! Default destructor:<br>
     */
    ~Screen();

    /*! Getter routines for members
     *  no setter routines are provided,
     *  these values are not expected to change once constructed */
    int get_Nx() { return Nx; }
    int get_Ny() { return Ny; }
    Vector get_Center()  { return Center; }
    Vector get_Normal()  { return Normal; }
    double get_dA()  { return dA; }
    Vector get_xVec()  { return xVec; }
    Vector get_yVec()  { return yVec; }
    
    /*! Get the position of one grid cell */
    Vector get_point(int ix, int iy);
    
    /*! Set the field trace data for one of the grid cells */
    void set_Trace(int ix, int iy, FieldTrace trace);
    
    /*! Total energy of radiation falling on the screen */
    double totalEnergy();
    
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
    
    /* array [Nx][Ny] of FieldTrace
     * outer (first) index is Nx
    */
    std::vector< std::vector<FieldTrace> > A;

};

