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

#include "screen.h"

Screen::Screen(
    int Nx_p, int Ny_p, int Nt_p,
    Vector xVec_p, Vector yVec_p,
    Vector center_p
    )
{
    Nx = Nx_p;
    Ny = Ny_p;
    xVec = xVec_p;
    yVec = yVec_p;
    Center = center_p;
    FieldTrace zero_trace(0.0,0.0,Nt_p);
    std::vector<FieldTrace> zero_column(Ny, zero_trace);
    A = std::vector< std::vector<FieldTrace> >(Nx,zero_column);
}

Screen::~Screen()
{
}

void Screen::set_Trace(int ix, int iy, FieldTrace trace)
{
    if (ix<0 || ix>=Nx) throw(Screen_IndexOutOfRange());
    if (iy<0 || iy>=Ny) throw(Screen_IndexOutOfRange());
    A[ix][iy] = trace;
}
