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

#include "vector.h"         // my own 3D vectors in space
#include <vector>           // the standard template library

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

    /* Poynting vector - energy flow density */
    Vector Poynting();
    
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

    /*! Division of the field by a real factor */
    ElMagField operator/ (double factor);

    /*! In-place multiplication of the field with a real factor */
    ElMagField& operator*= (double factor);

    /*! in-place division of the field by a real factor */
    ElMagField& operator/= (double factor);

private:

    Vector vecE;
    Vector vecB;

};

/*  Assignments of field traces may go wrong -> throws an exception */
class FieldTrace_SizeMismatch { };
class FieldTrace_IndexOutOfRange { };
class FieldTrace_Zero { };

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

    /*! Assignment operator which copies all data of another field trace.
     *  An exception is thrown if the sizes don't match.
     */
    FieldTrace & operator=(const FieldTrace &t);

    /*! Set one of the entries of a field trace */
    void set(int index, ElMagField f);

    /*! Set all entries of a field trace from a buffer.
     *  The number of elements in the buffer needs to be provided for checking.
     */
    void set(ElMagField *buffer, int Nb);

    /*! Multiplication of the field with a real factor */
    FieldTrace operator* (double factor);
    
    /*! Division of the field by a real factor */
    FieldTrace operator/ (double factor);
    
    /*! Sum of the fields of two traces.
     *  An exception is thrown in case of mismatch of the time or size definitions.
     */
    FieldTrace operator+ (FieldTrace other);

    /*! Accumulating sum of the fields of two traces.
     *  An exception is thrown in case of mismatch of the time or size definitions.
     */
    FieldTrace& operator+= (FieldTrace other);

    /*! Difference of the fields of two traces.
     *  An exception is thrown in case of mismatch of the time or size definitions.
     */
    FieldTrace operator- (FieldTrace other);

    /*! Get one of the sampling definitions
     *  These are not expected to be changed, so, no setter routines are provoded */
    double get_t0() { return t0; };
    double get_dt() { return dt; };
    int get_N() { return N; };
    
    /*! Get the time of one point of the trace */
    double get_time(int index);
    
    /*! Get the field at one point of the trace */
    ElMagField get_field(int index);
    
    /*! Get the field at one point in time.
        This interpolates linearly between samples.
        Outside the covered trace length zero is returned.
     */
    ElMagField get_field(double time);
    
    /*! Poynting vector - time-integrated energy flow density */
    Vector Poynting();
    
    /*! compute time derivative of fields */
    FieldTrace derivative();
    
    /*! compute a retarded trace
     *  @param delta_t - the amount of retardation in seconds
     *  @param t0_p - start time of the new trace in seconds
     *  @param N_p - length of the new trace
     *  the time step remains the same
     */
    FieldTrace retarded(double delta_t, double t0_p, int N_p);

private:

    int N;
    double t0;
    double dt;
    std::vector<ElMagField> trace;

};

