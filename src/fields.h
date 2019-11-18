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

/*!
 * \class ElMagField
 * \brief Type for electromagnetic fields.
 * @author Ulf Lehnert
 * @date 18.11.2019
 * 
 * This class holds the two Vectors E [V/m] and B [T] of the electromagnetic field.
 * Memory is 6 Doubles for the two Vectors.
 */
class ElMagField
{

public:

    /*! Default constructor:<br>
     * initializes both components to zero
     */
    ElMagField();

    /*! Initializing constructor: */
    ElMagField(Vector E, Vector B);

    /*! Set all components to zero */
    void Zero();

    /*! Electric field report */
    Vector E();

    /*! Magnetic field report */
    Vector B();

    /*! Sum of two fields */
    ElMagField operator+ (ElMagField other);

    /*! Accumulating sum of two fields */
    ElMagField& operator+= (ElMagField other);

    /*! Difference of two fields */
    ElMagField operator- (ElMagField other);

    /*! In-place difference of two fields */
    ElMagField& operator-= (ElMagField other);

    /*! Multiplication of the field with a real factor */
    ElMagField operator* (double factor);

    /*! In-place multiplication of the field with a real factor */
    ElMagField& operator*= (double factor);

private:

    Vector vecE;
    Vector vecB;

};

/*!
 * \class FieldTrace
 * \brief Type for time traces of electromagnetic fields.
 * @author Ulf Lehnert
 * @date 18.11.2019
 * 
 * This class holds a number of the electromagnetic field samples
 * covering a certain time window. Time of the first sample is t0.
 * All other samples are spaced by dt from each other.
 */
class FieldTrace
{

public:

    /*! Default constructor:<br>
     *  Allocate memory for N samples
     *  and initialize all components to zero.
     */
    FieldTrace(double t0_p, double dt_p, int N_p);

    /*! Default destructor:<br>
     *  Free memory
     */
    ~FieldTrace();

    double get_t0() { return t0; };
    double get_dt() { return dt; };
    int get_N() { return N; };
    
private:

    int N;
    double t0;
    double dt;
    ElMagField *trace;

};

