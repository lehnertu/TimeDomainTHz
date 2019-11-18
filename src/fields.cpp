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

#include "fields.h"

// ********** ElMagField **********

ElMagField::ElMagField()
{
    Zero();
}

ElMagField::ElMagField(Vector E, Vector B)
{
    vecE=E;
    vecB=B;
}

void ElMagField::Zero()
{
    vecE=Vector(0.0,0.0,0.0);
    vecB=Vector(0.0,0.0,0.0);
}

Vector ElMagField::E() { return vecE; }

Vector ElMagField::B() { return vecB; }

ElMagField ElMagField::operator+ (ElMagField other)
{
    ElMagField temp;
    temp.vecE = vecE + other.vecE;
    temp.vecB = vecB + other.vecB;
    return (temp);
}

ElMagField& ElMagField::operator+= (ElMagField other)
{
    vecE += other.vecE;
    vecB += other.vecB;
    return (*this);
}

ElMagField ElMagField::operator- (ElMagField other)
{
    ElMagField temp;
    temp.vecE = vecE - other.vecE;
    temp.vecB = vecB - other.vecB;
    return (temp);
}

ElMagField& ElMagField::operator-= (ElMagField other)
{
    vecE -= other.vecE;
    vecB -= other.vecB;
    return (*this);
}

ElMagField ElMagField::operator* (double factor)
{
    ElMagField temp;
    temp.vecE = vecE*factor;
    temp.vecB = vecB*factor;
    return (temp);
}

ElMagField& ElMagField::operator*= (double factor)
{
    vecE *= factor;
    vecB *= factor;
    return (*this);
}

FieldTrace::FieldTrace(double t0_p, double dt_p, int N_p)
{
    t0 = t0_p;
    dt = dt_p;
    N = N_p;
    trace = new ElMagField[N];
}

FieldTrace::~FieldTrace()
{
    delete [] trace;
}


