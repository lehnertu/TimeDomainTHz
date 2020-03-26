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
#include "global.h"
#include <exception>
#include <math.h>
#include <iostream>

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

Vector ElMagField::Poynting() { return cross(vecE, vecB) / MuNull; }

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

ElMagField ElMagField::operator/ (double factor)
{
    ElMagField temp;
    temp.vecE = vecE/factor;
    temp.vecB = vecB/factor;
    return (temp);
}

ElMagField& ElMagField::operator*= (double factor)
{
    vecE *= factor;
    vecB *= factor;
    return (*this);
}

ElMagField& ElMagField::operator/= (double factor)
{
    vecE /= factor;
    vecB /= factor;
    return (*this);
}

// ========================================================
//    FieldTrace methods
// ========================================================

FieldTrace::FieldTrace(double t0_p, double dt_p, int N_p)
{
    t0 = t0_p;
    dt = dt_p;
    N = N_p;
    trace = std::vector<ElMagField>(N,ElMagField(Vector(0.0,0.0,0.0),Vector(0.0,0.0,0.0)));
}

FieldTrace::FieldTrace(FieldTrace *original)
{
    t0 = original->t0;
    dt = original->dt;
    N = original->N;
    trace = original->trace;
}

void FieldTrace::zero()
{
    for (int i=0; i<N; i++)
        trace[i] = ElMagField(Vector(0.0,0.0,0.0),Vector(0.0,0.0,0.0));
}

FieldTrace::~FieldTrace()
{
    // delete [] trace;
}

FieldTrace & FieldTrace::operator=(const FieldTrace &t)
{
    // check for "self assignment" and do nothing in that case
    if (this != &t)
    {
        // check for a size mismatch and throw an exception in that case
        if (N != t.N)
        {
            throw(FieldTrace_SizeMismatch());
        }
        // copy all data from the other field trace
        else
        {
            t0 = t.t0;
            dt = t.dt;
            // for (int i=0; i<N; i++) trace[i]=t.trace[i];
            trace = t.trace;
        }
    }
    return *this;    
}

FieldTrace FieldTrace::operator* (double factor)
{
    FieldTrace temp = *this;
    for (int i=0; i<temp.N; i++)
        temp.trace[i] *= factor;
    return(temp);
}

FieldTrace FieldTrace::operator/ (double factor)
{
    FieldTrace temp = *this;
    for (int i=0; i<temp.N; i++)
        temp.trace[i] /= factor;
    return(temp);
}

FieldTrace FieldTrace::operator+ (FieldTrace other)
{
    if (N != other.N) throw(FieldTrace_SizeMismatch());
    if (t0 != other.t0) throw(FieldTrace_TimeMismatch());
    if (dt != other.dt) throw(FieldTrace_TimeMismatch());
    FieldTrace temp = *this;
    for (int i=0; i<N; i++)
        temp.trace[i] = trace[i] + other.trace[i];
    return (temp);
}

FieldTrace& FieldTrace::operator+= (FieldTrace other)
{
    if (N != other.N) throw(FieldTrace_SizeMismatch());
    if (t0 != other.t0) throw(FieldTrace_TimeMismatch());
    if (dt != other.dt) throw(FieldTrace_TimeMismatch());
    for (int i=0; i<N; i++)
        trace[i] += other.trace[i];
    return (*this);
}

FieldTrace FieldTrace::operator- (FieldTrace other)
{
    if (N != other.N) throw(FieldTrace_SizeMismatch());
    if (t0 != other.t0) throw(FieldTrace_TimeMismatch());
    if (dt != other.dt) throw(FieldTrace_TimeMismatch());
    FieldTrace temp = *this;
    for (int i=0; i<N; i++)
        temp.trace[i] = trace[i] - other.trace[i];
    return (temp);
}

void FieldTrace::set(int index, ElMagField f)
{
    if (index<0 || index>=N) throw(FieldTrace_IndexOutOfRange());
    trace[index] = f;
}

void FieldTrace::set(ElMagField *buffer, int Nb)
{
    if (Nb != N) throw(FieldTrace_SizeMismatch());
    ElMagField *buf = buffer;
    for (int it=0; it<N; it++) set(it, *buf++);
}

double FieldTrace::get_time(int index)
{
    if (index<0 || index>=N) throw(FieldTrace_IndexOutOfRange());
    return t0+index*dt;
}

ElMagField FieldTrace::get_field(int index)
{
    if (index<0 || index>=N) throw(FieldTrace_IndexOutOfRange());
    return trace[index];
}

ElMagField FieldTrace::get_field(double time)
{
    double ref = (time-t0)/dt;
    int index = (int)floor(ref);
    double frac = ref-(double)index;
    ElMagField f1(Vector(0.0,0.0,0.0),Vector(0.0,0.0,0.0));
    if ((index>=0) && (index<N)) f1=trace[index];
    ElMagField f2(Vector(0.0,0.0,0.0),Vector(0.0,0.0,0.0));
    if ((index+1>=0) && (index+1<N)) f2=trace[index+1];
    return(f1*(1.0-frac)+f2*frac);
}

Vector FieldTrace::Poynting()
{
    Vector sum(0.0,0.0,0.0);
    for (int i=0; i<N; i++)
        sum += trace[i].Poynting();
    return(sum*dt);
}

FieldTrace FieldTrace::derivative()
{
    FieldTrace temp = *this;
    temp.trace[0] = (trace[1]-trace[0])/dt;
    for (int i=1; i<N-1; i++)
        temp.trace[i] = (trace[i+1]-trace[i-1])/(2.0*dt);
    temp.trace[N-1] = (trace[N-1]-trace[N-2])/dt;
    return(temp);
}

void FieldTrace::retard(double delta_t, FieldTrace *target)
{
    for (int i=0; i<target->get_N(); i++)
    {
        double t = target->get_time(i) - delta_t;
        target->trace[i] = get_field(t);
    };
    if (target->Poynting().norm()==0.0)
    {
        if (DEBUGLEVEL>=2)
        { 
            std::cout << "retard from t0=" << t0 << " N=" << N << " P=" << Poynting().norm() << std::endl;
            std::cout << "      by delta=" << delta_t << "s" << std::endl;
            std::cout << "retard   to t0=" << target->get_t0() << " N=" << target->get_N();
            std::cout << " P=" << target->Poynting().norm() << std::endl;
        }
        throw(FieldTrace_Zero());
    };
}

void FieldTrace::transform(Vector ex, Vector ey, Vector ez)
{
    for (int i=0; i<N; i++)
    {
        Vector E = trace[i].E();
        Vector B = trace[i].B();
        Vector E_new = Vector(dot(E,ex),dot(E,ey),dot(E,ez));
        Vector B_new = Vector(dot(B,ex),dot(B,ey),dot(B,ez));
        trace[i] = ElMagField(E_new, B_new);
    }
}

