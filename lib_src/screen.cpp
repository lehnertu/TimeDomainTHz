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

#include "global.h"
#include "screen.h"
#include "hdf5.h"

Screen::Screen(
    int Ncp_p, int Np_p, int Nt_p,
    double dt_p
    )
{
    Ncp = Ncp_p;
    Np = Np_p;
    Nt = Nt_p;
    dt = dt_p;
    
    triangle_points = std::vector<Vector>(Ncp);
    field_points = std::vector<Vector>(Np);
    t0 = std::vector<double>(Nt,0.0);
    FieldTrace zero_trace(0.0,dt,Nt);
    A = std::vector<FieldTrace>(Np,zero_trace);
}

Screen::Screen(std::string filename)
{
    herr_t status;
    hid_t dataset, attr;
    cout << "reading HDF5 file " << filename << endl;
    hid_t file = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file<0) throw(Screen_FileReadError());
    
    // read the mesh corner points
    dataset = H5Dopen2(file, "MeshCornerPoints", H5P_DEFAULT);
    if (dataset<0) throw(Screen_FileReadError());
    attr = H5Aopen_name(dataset, "Ncp");
    if (attr<0) throw(Screen_FileReadError());
    status =  H5Aread(attr, H5T_NATIVE_INT, &Ncp);
    if (status<0) throw(Screen_FileReadError());
    Vector *pos = new Vector[Ncp];
    status = H5Dread (dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        pos);
    if (status<0) throw(Screen_FileReadError());
    status = H5Dclose(dataset);
    if (status<0) throw(Screen_FileReadError());
    triangle_points = std::vector<Vector>(Ncp);
    for (int icp=0; icp<Ncp; icp++)
        triangle_points[icp] = pos[icp];
    delete pos;
    
    // read the mesh triangles
    
    // read the position dataset
    dataset = H5Dopen2(file, "ObservationPosition", H5P_DEFAULT);
    if (dataset<0) throw(Screen_FileReadError());
    attr = H5Aopen_name(dataset, "Np");
    if (attr<0) throw(Screen_FileReadError());
    status =  H5Aread(attr, H5T_NATIVE_INT, &Np);
    if (status<0) throw(Screen_FileReadError());
    // create and fill a buffer
    pos = new Vector[Np];
    status = H5Dread (dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        pos);
    if (status<0) throw(Screen_FileReadError());
    status = H5Dclose(dataset);
    if (status<0) throw(Screen_FileReadError());
    field_points = std::vector<Vector>(Np);
    for (int ip=0; ip<Np; ip++)
        field_points[ip] = pos[ip];
    delete pos;

    // read the trace timing dataset
    // it is not an error if there are no data available
    dataset = H5Dopen2(file, "ObservationTime", H5P_DEFAULT);
    if (dataset<0)
    {
        // no data available
        Nt = 0;
        dt = 0.0;
        t0 = std::vector<double>(Np,0.0);
    }
    else
    {
        attr = H5Aopen_name(dataset, "Nt");
        if (attr<0) throw(Screen_FileReadError());
        status =  H5Aread(attr, H5T_NATIVE_INT, &Nt);
        if (status<0) throw(Screen_FileReadError());
        attr = H5Aopen_name(dataset, "dt");
        if (attr<0) throw(Screen_FileReadError());
        status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &dt);
        if (status<0) throw(Screen_FileReadError());
        double *t0_buf = new double[Np];
        status = H5Dread (dataset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            H5S_ALL,
            H5P_DEFAULT,			// data transfer properties
            t0_buf);
        if (status<0) throw(Screen_FileReadError());
        status = H5Dclose(dataset);
        if (status<0) throw(Screen_FileReadError());
        t0 = std::vector<double>(Np);
        for (int ip=0; ip<Np; ip++)
            t0[ip] = t0_buf[ip];
        delete t0_buf;
    }

    // read the field dataset
    // it is not an error if there are no data available
    dataset = H5Dopen2(file, "ElMagField", H5P_DEFAULT);
    if (dataset<0)
    {
        // no data available
        FieldTrace zero_trace(0.0,dt,Nt);
        A = std::vector<FieldTrace>(Np,zero_trace);
    }
    else
    {
        ElMagField *field = new ElMagField[Np*Nt];
        status = H5Dread (dataset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            H5S_ALL,
            H5P_DEFAULT,			// data transfer properties
            field);
        // transfer the data into the internal data structures
        FieldTrace trace(0.0, dt, Nt);
        A = std::vector<FieldTrace>(Np,trace);
        ElMagField *buf = field;
        for (int ip=0; ip<Np; ip++)
        {   
            trace.set_t0(t0[ip]);
            trace.set(buf,Nt);
            buf += Nt;
            A[ip] = trace;
        }
        delete field;
    }
        
    status = H5Fclose (file);
    if (status<0) throw(Screen_FileReadError());
}

Screen::~Screen()
{
}

void Screen::writeReport(std::ostream *st)
{
    *st << "TimeDomainTHz - Screen" << std::endl << std::endl;
    *st << "Np=" << Np << "  Nt=" << get_Nt() << std::endl;
    /*
    st->precision(4);
    *st << "center = (" << Center.x << ", " << Center.y << ", " << Center.z << ") m" << std::endl;
    *st << "e_x = (" << xVec.x*1.0e3 << ", " << xVec.y*1.0e3 << ", " << xVec.z*1.0e3 << ") mm" << std::endl;
    *st << "e_y = (" << yVec.x*1.0e3 << ", " << yVec.y*1.0e3 << ", " << yVec.z*1.0e3 << ") mm" << std::endl;
    *st << "n   = (" << Normal.x << ", " << Normal.y << ", " << Normal.z << ")" << std::endl;
    st->precision(6);
    *st << "dA  = " << dA*1.0e6 << " mm²" << std::endl;
    *st << "t0  = " << A[Nx/2][Ny/2].get_t0()*1.0e9 << " ns" << std::endl;
    *st << "dt  = " << get_dt()*1.0e9 << " ns" << std::endl << std::endl;
    */
    int peak_index=-1;
    double peak=0.0;
    Vector Sp;
    for (int ip=0; ip<Np; ip++)
    {   
        Vector p = A[ip].Poynting();
        if (p.norm()>peak)
        {
            Sp = p;
            peak = p.norm();
            peak_index = ip;
        }
    }
    st->precision(6);
    // *st << "Energy incident on screen = " << totalEnergy()*1.0e6 << " µJ" << std::endl << std::endl;
    Vector p = get_point(peak_index);
    *st << "pos["  << peak_index << "] =   (" << p.x << ", " << p.y << ", " << p.z << ") m   " << std::endl;
    *st << "peak energy density =" << peak << " J/m²" << std::endl;
    *st << "Poynting vector S = (" << Sp.x << ", " << Sp.y << ", " << Sp.z << ") J/m²" << std::endl;
    *st << "t0=" << A[peak_index].get_t0()*1e9 << " ns  " << std::endl;
    *st << std::endl;
}

void Screen::writeTraceReport(std::ostream *st, int ip)
{
    *st << "trace[" << ip << "]  ";
    Vector p = get_point(ip);
    *st << "pos=(" << p.x << ", " << p.y << ", " << p.z << ") m   " << std::endl;
    *st << "t0=" << A[ip].get_t0()*1e9 << " ns  ";
    *st << "P=" << A[ip].Poynting().norm() << " J/m²" << std::endl;
}

