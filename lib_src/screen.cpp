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
    int Nx_p, int Ny_p, int Nt_p,
    Vector xVec_p, Vector yVec_p,
    Vector center_p
    )
{
    Nx = Nx_p;
    Ny = Ny_p;
    xVec = xVec_p;
    yVec = yVec_p;
    Normal = cross(xVec,yVec);
    dA = Normal.norm();
    if (dot(xVec,yVec)>1e-6*dA) throw(Screen_NonOrthogonal());
    Normal.normalize();
    Center = center_p;
    FieldTrace zero_trace(0.0,0.0,Nt_p);
    std::vector<FieldTrace> zero_column(Ny, zero_trace);
    A = std::vector< std::vector<FieldTrace> >(Nx,zero_column);
    dn_A = std::vector< std::vector<FieldTrace> >(Nx,zero_column);
    dt_A = std::vector< std::vector<FieldTrace> >(Nx,zero_column);
}

Screen::Screen(std::string filename)
{
    herr_t status;
    hid_t attr;
    cout << "reading HDF5 file " << filename << endl;
    hid_t file = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file<0) throw(Screen_FileReadError());
    
    // open the position dataset
    hid_t pos_dataset = H5Dopen2(file, "ObservationPosition", H5P_DEFAULT);
    if (pos_dataset<0) throw(Screen_FileReadError());
    attr = H5Aopen_name(pos_dataset, "Nx");
    if (attr<0) throw(Screen_FileReadError());
    status =  H5Aread(attr, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(Screen_FileReadError());
    attr = H5Aopen_name(pos_dataset, "Ny");
    if (attr<0) throw(Screen_FileReadError());
    status =  H5Aread(attr, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(Screen_FileReadError());
    // create and fill a buffer
    Vector *pos = new Vector[Nx*Ny];
    status = H5Dread (pos_dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        pos);
    if (status<0) throw(Screen_FileReadError());
    status = H5Dclose(pos_dataset);
    if (status<0) throw(Screen_FileReadError());
    // analyze the data
    Vector pos00 = pos[0];
    Vector pos01 = pos[Ny-1];
    Vector pos10 = pos[(Nx-1)*Ny];
    Vector pos11 = pos[(Nx-1)*Ny+Ny-1];
    xVec = (pos10-pos00)/(double)(Nx-1);
    yVec = (pos01-pos00)/(double)(Ny-1);
    Normal = cross(xVec,yVec);
    dA = Normal.norm();
    if (dot(xVec,yVec)>1e-6*dA) throw(Screen_NonOrthogonal());
    Normal.normalize();
    Center = (pos11+pos00)*0.5;
    delete pos;

    // open the trace timing dataset
    hid_t time_dataset = H5Dopen2(file, "ObservationTime", H5P_DEFAULT);
    if (time_dataset<0) throw(Screen_FileReadError());
    attr = H5Aopen_name(time_dataset, "Nt");
    if (attr<0) throw(Screen_FileReadError());
    int Nt;
    status =  H5Aread(attr, H5T_NATIVE_INT, &Nt);
    if (status<0) throw(Screen_FileReadError());
    attr = H5Aopen_name(time_dataset, "dt");
    if (attr<0) throw(Screen_FileReadError());
    double dt;
    status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &dt);
    if (status<0) throw(Screen_FileReadError());
    // create and fill a buffer
    double *t0_buf = new double[Nx*Ny];
    status = H5Dread (time_dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        t0_buf);
    if (status<0) throw(Screen_FileReadError());
    status = H5Dclose(time_dataset);
    if (status<0) throw(Screen_FileReadError());

    // open the field dataset
    hid_t field_dataset = H5Dopen2(file, "ElMagField", H5P_DEFAULT);
    if (field_dataset<0) throw(Screen_FileReadError());
    // create and fill a buffer
    ElMagField *field = new ElMagField[Nx*Ny*Nt];
    status = H5Dread (field_dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        field);
    // transfer the data into the internal data structures
    FieldTrace trace(0.0, dt, Nt);
    std::vector<FieldTrace> zero_column(Ny, trace);
    A = std::vector< std::vector<FieldTrace> >(Nx,zero_column);
    ElMagField *buf = field;
    double *pt0 = t0_buf;
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {   
            trace.set_t0(*pt0++);
            trace.set(buf,Nt);
            buf += Nt;
            A[ix][iy] = trace;
        }
    delete t0_buf;
    delete field;
    
    status = H5Fclose (file);
    if (status<0) throw(Screen_FileReadError());
    dn_A = std::vector< std::vector<FieldTrace> >(Nx,zero_column);
    dt_A = std::vector< std::vector<FieldTrace> >(Nx,zero_column);
}

Screen::~Screen()
{
}

Vector Screen::get_point(int ix, int iy)
{
    return Center + xVec*((double)ix-0.5*((double)Nx-1.0)) + yVec*((double)iy-0.5*((double)Ny-1.0));
}

FieldTrace Screen::get_trace(int ix, int iy)
{
    if (ix<0 || ix>=Nx) throw(Screen_IndexOutOfRange());
    if (iy<0 || iy>=Ny) throw(Screen_IndexOutOfRange());
    return A[ix][iy];
}

void Screen::set_trace(int ix, int iy, FieldTrace trace)
{
    if (ix<0 || ix>=Nx) throw(Screen_IndexOutOfRange());
    if (iy<0 || iy>=Ny) throw(Screen_IndexOutOfRange());
    A[ix][iy] = trace;
}

double Screen::totalEnergy()
{
    double sum = 0.0;
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
            sum += -dot(A[ix][iy].Poynting(),Normal);
    return(sum*dA);
}

void Screen::writeReport(std::ostream *st)
{
    *st << "TimeDomainTHz - Screen" << std::endl << std::endl;
    *st << "Nx=" << Nx << "  Ny=" << Ny << "  Nt=" << get_Nt() << std::endl;
    st->precision(4);
    *st << "center = (" << Center.x << ", " << Center.y << ", " << Center.z << ") m" << std::endl;
    *st << "e_x = (" << xVec.x*1.0e3 << ", " << xVec.y*1.0e3 << ", " << xVec.z*1.0e3 << ") mm" << std::endl;
    *st << "e_y = (" << yVec.x*1.0e3 << ", " << yVec.y*1.0e3 << ", " << yVec.z*1.0e3 << ") mm" << std::endl;
    *st << "n   = (" << Normal.x << ", " << Normal.y << ", " << Normal.z << ")" << std::endl;
    st->precision(6);
    *st << "dA  = " << dA*1.0e6 << " mm²" << std::endl;
    *st << "t0  = " << A[Nx/2][Ny/2].get_t0()*1.0e9 << " ns" << std::endl;
    *st << "dt  = " << get_dt()*1.0e9 << " ns" << std::endl << std::endl;
    double pp=0.0;
    Vector Sp;
    int pix=-1;
    int piy=-1;
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {   
            Vector p=A[ix][iy].Poynting();
            if (p.norm()>pp)
            {
                Sp = p;
                pp = p.norm();
                pix=ix;
                piy=iy;
            }
        }
    st->precision(6);
    *st << "peak energy density (" << pix << ", " << piy << ") = " << pp << " J/m²" << std::endl;
    *st << "Energy incident on screen = " << totalEnergy()*1.0e6 << " µJ" << std::endl << std::endl;
}

void Screen::writeTraceReport(std::ostream *st, int ix, int iy)
{
    *st << "trace (" << ix << ", " << iy << ")  ";
    Vector p = get_point(ix,iy);
    *st << "pos=(" << p.x << ", " << p.y << ", " << p.z << ") m   ";
    *st << "t0=" << A[ix][iy].get_t0()*1e9 << " ns  ";
    *st << "P=" << A[ix][iy].Poynting().norm() << " J/m²" << std::endl;
}

void Screen::bufferArray(double *buffer)
{
    double* bp = buffer;
    int Nt = A[0][0].get_N();
    FieldTrace trace(0.0,0.0,Nt);
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            trace = A[ix][iy];
            // for testing purposes this line can be replaced like
            // trace = dx_A[ix][iy];
            // then the computed derivative is written to an output file instead of the field
            for (int it=0; it<Nt; it++)
            {
                ElMagField field = trace.get_field(it);
                *bp++ = field.E().x;
                *bp++ = field.E().y;
                *bp++ = field.E().z;
                *bp++ = field.B().x;
                *bp++ = field.B().y;
                *bp++ = field.B().z;
            };
        }
}

FieldTrace Screen::dx_A(int ix, int iy)
{
    if (ix<0 || ix>=Nx) throw(Screen_IndexOutOfRange());
    if (iy<0 || iy>=Ny) throw(Screen_IndexOutOfRange());
    double dX = xVec.norm();
    FieldTrace center = A[ix][iy];
    // the neighbouring traces are resampled to the timing defined by the center trace
    FieldTrace right = center;
    FieldTrace left = center;
    // the return trace
    FieldTrace trace = center;
    if (ix==0)
    {
        A[1][iy].retard(0.0, &right);
        trace = (right-center) / dX;
    }
    else 
    {
        if (ix==Nx-1)
        {
            A[Nx-2][iy].retard(0.0, &left);
            trace = (center-left) / dX;
        }
        else
        {
            A[ix+1][iy].retard(0.0, &right);
            A[ix-1][iy].retard(0.0, &left);
            trace = (right-left) / (dX*2.0);
        };
    }
    return trace;
}

FieldTrace Screen::dy_A(int ix, int iy)
{
    if (ix<0 || ix>=Nx) throw(Screen_IndexOutOfRange());
    if (iy<0 || iy>=Ny) throw(Screen_IndexOutOfRange());
    double dY = yVec.norm();
    FieldTrace center = A[ix][iy];
    // the neighbouring traces are resampled to the timing defined by the center trace
    FieldTrace right = center;
    FieldTrace left = center;
    // the return trace
    FieldTrace trace = center;
    if (iy==0)
    {
        A[ix][1].retard(0.0, &right);
        trace = (right-center) / dY;
    }
    else 
    {
        if (iy==Ny-1)
        {
            A[ix][Ny-2].retard(0.0, &left);
            trace = (center-left) / dY;
        }
        else
        {
            A[ix][iy+1].retard(0.0, &right);
            A[ix][iy-1].retard(0.0, &left);
            trace = (right-left) / (dY*2.0);
        };
    }
    return trace;
}

void Screen::computeDerivatives()
{
    int Nt = A[0][0].get_N();
    FieldTrace trace(0.0,0.0,Nt);
    // All computations are done in the screen-local reference frame.
    // Is is assumed that (xVec, yVec, Normal) form a right-handed orthogonal reference frame.
    // Components Ex, Bx are those in the xVec direction (unit vector ex)
    // and all other components accordingly
    Vector ex = xVec;
    ex.normalize();
    Vector ey = yVec;
    ey.normalize();
    Vector ez = Normal;
    // declare some pointers needed to copy field components
    double *pEx, *pEy, *pEz, *pBx, *pBy, *pBz;
    // compute the time-derivatives for all traces
    // this derivative is stored in the global reference frame
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
            dt_A[ix][iy] = A[ix][iy].derivative();
    // compute the screen-local components of the time derivatives
    double *dt_Ex = new double[Nx*Ny*Nt];
    double *dt_Ey = new double[Nx*Ny*Nt];
    double *dt_Bx = new double[Nx*Ny*Nt];
    double *dt_By = new double[Nx*Ny*Nt];
    pEx = dt_Ex;
    pEy = dt_Ey;
    pBx = dt_Bx;
    pBy = dt_By;
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            // here the components of the field are in global coordinates
            trace = dt_A[ix][iy];
            // translate to screen-local components
            for (int it=0; it<Nt; it++)
            {
                ElMagField f = trace.get_field(it);
                *pEx++ = dot(f.E(),ex);
                *pEy++ = dot(f.E(),ey);
                *pBx++ = dot(f.B(),ex);
                *pBy++ = dot(f.B(),ey);
            }
        };
    // Compute the derivatives along the xVec axis
    double *dx_Ex = new double[Nx*Ny*Nt];
    double *dx_Ez = new double[Nx*Ny*Nt];
    double *dx_Bx = new double[Nx*Ny*Nt];
    double *dx_Bz = new double[Nx*Ny*Nt];
    pEx = dx_Ex;
    pEz = dx_Ez;
    pBx = dx_Bx;
    pBz = dx_Bz;
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            // here the components of the field are in global coordinates
            trace = dx_A(ix,iy);
            // translate to screen-local components
            for (int it=0; it<Nt; it++)
            {
                ElMagField f = trace.get_field(it);
                *pEx++ = dot(f.E(),ex);
                *pEz++ = dot(f.E(),ez);
                *pBx++ = dot(f.B(),ex);
                *pBz++ = dot(f.B(),ez);
            }
        };
    // Compute the derivatives along the yVec axis
    double *dy_Ey = new double[Nx*Ny*Nt];
    double *dy_Ez = new double[Nx*Ny*Nt];
    double *dy_By = new double[Nx*Ny*Nt];
    double *dy_Bz = new double[Nx*Ny*Nt];
    pEy = dy_Ey;
    pEz = dy_Ez;
    pBy = dy_By;
    pBz = dy_Bz;
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            // here the components of the field are in global coordinates
            trace = dy_A(ix,iy);
            // translate to screen-local components
            for (int it=0; it<Nt; it++)
            {
                ElMagField f = trace.get_field(it);
                *pEy++ = dot(f.E(),ey);
                *pEz++ = dot(f.E(),ez);
                *pBy++ = dot(f.B(),ey);
                *pBz++ = dot(f.B(),ez);
            }
        };
    // Compute the ez derivatives using the Maxwell equations in screen-local coordinates.
    // This is the drivative along the Normal direction (in the screen-local frame).
    double *dz_Ex = new double[Nx*Ny*Nt];
    double *dz_Ey = new double[Nx*Ny*Nt];
    double *dz_Ez = new double[Nx*Ny*Nt];
    double *dz_Bx = new double[Nx*Ny*Nt];
    double *dz_By = new double[Nx*Ny*Nt];
    double *dz_Bz = new double[Nx*Ny*Nt];
    double cSquared = SpeedOfLight*SpeedOfLight;
    for (int i=0; i<Nx*Ny*Nt; i++)
    {
        dz_Ex[i] = dx_Ez[i] - dt_By[i];
        dz_Ey[i] = dy_Ez[i] + dt_Bx[i];
        dz_Ez[i] = -dx_Ex[i] - dy_Ey[i];
        dz_Bx[i] = dx_Bz[i] + dt_Ey[i]/cSquared;
        dz_By[i] = dy_Bz[i] - dt_Ex[i]/cSquared;
        dz_Bz[i] = -dx_Bx[i] - dy_By[i];
    };
    // Transform the normal derivatives back to global coordinates
    pEx = dz_Ex;
    pEy = dz_Ey;
    pEz = dz_Ez;
    pBx = dz_Bx;
    pBy = dz_By;
    pBz = dz_Bz;
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            // this is necessary to set the correct timing
            trace = A[ix][iy];
            // now set the field data
            for (int it=0; it<Nt; it++)
            {
                Vector E = ex * (*pEx++) + ey * (*pEy++) + ez * (*pEz++);
                Vector B = ex * (*pBx++) + ey * (*pBy++) + ez * (*pBz++);
                ElMagField f = ElMagField(E,B);
                trace.set(it,f);
            }
            dn_A[ix][iy] = trace;
        };
    delete dx_Ex;
    delete dx_Ez;
    delete dx_Bx;
    delete dx_Bz;
    delete dy_Ey;
    delete dy_Ez;
    delete dy_By;
    delete dy_Bz;
    delete dt_Ex;
    delete dt_Ey;
    delete dt_Bx;
    delete dt_By;
    delete dz_Ex;
    delete dz_Ey;
    delete dz_Ez;
    delete dz_Bx;
    delete dz_By;
    delete dz_Bz;
}

void Screen::propagate_to(Vector target_pos, FieldTrace *target_trace)
{
    // this is the result sum
    FieldTrace trace(target_trace);
    trace.zero();
    // this is the source trace
    FieldTrace t1(A[0][0]);
    // this is the component to be added to the result trace
    FieldTrace t2(target_trace);
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {   
            Vector source = get_point(ix,iy);
            Vector RVec = target_pos - source;
            double R = RVec.norm();
            double R2 = R*R;
            double R3 = R2*R;
            try
            {
                t1 = A[ix][iy];
                t1.retard(R/SpeedOfLight, &t2);
                trace += t2 * (-dot(RVec,Normal)/R3);
            }
            catch(FieldTrace_Zero exc)
            {
                std::cout << "FieldTrace first term zero exception at ix=" << ix << "  iy=" << iy << std::endl;
                throw(Screen_Zero_exception(ix,iy));
            }
            try
            {
                t1 = dt_A[ix][iy];
                t1.retard(R/SpeedOfLight, &t2);
                trace += t2 * (-dot(RVec,Normal)/(R2*SpeedOfLight));
            }
            catch(FieldTrace_Zero exc)
            {
                std::cout << "FieldTrace second term zero exception at ix=" << ix << "  iy=" << iy << std::endl;
                throw(Screen_Zero_exception(ix,iy));
            }
            try
            {
                t1 = dn_A[ix][iy];
                t1.retard(R/SpeedOfLight, &t2);
                // TODO: why is this term positive - should be negative
                trace += t2 * (1.0/R);
            }
            catch(FieldTrace_Zero exc)
            {
                std::cout << "FieldTrace third term zero exception at ix=" << ix << "  iy=" << iy << std::endl;
                throw(Screen_Zero_exception(ix,iy));
            }
        }    
    *target_trace = trace * dA/(4.0*Pi);
}

void Screen::writeFieldHDF5(std::string filename)
{
    herr_t status;
    hid_t atts, att;
    hid_t dataspace, dcpl, dataset;
    
    cout << "writing HDF5 file " << filename << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file<0) throw(Screen_FileWriteError());
    
    // Create dataspace for the observation positions.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t pdims[3];
    pdims[0] = Nx;
    pdims[1] = Ny;
    pdims[2] = 3;
    dataspace = H5Screate_simple (3, pdims, NULL);
    if (dataspace<0) throw(Screen_FileWriteError());
    // buffer the data
    double *buffer = new double[Nx*Ny*3];
    double *bp = buffer;
    for (int ix=0; ix<Nx; ix++)
	    for (int iy=0; iy<Ny; iy++)
    	{
	        Vector pos = get_point(ix, iy);
	        *bp++ = pos.x;
	        *bp++ = pos.y;
	        *bp++ = pos.z;
	    };
    // Create the dataset creation property list
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(Screen_FileWriteError());
    // Create the dataset.
    dataset = H5Dcreate(file,
        "ObservationPosition",		// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        dataspace, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dataset<0) throw(Screen_FileWriteError());
    // Write the data to the dataset
    status = H5Dwrite (dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        dataspace,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(Screen_FileWriteError());
    // attach scalar attributes
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "Nx", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    status = H5Awrite(att, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "Ny", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    status = H5Awrite(att, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Sclose (atts);
    if (status<0) throw(Screen_FileWriteError());
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Dclose (dataset);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Sclose (dataspace);
    if (status<0) throw(Screen_FileWriteError());
    delete buffer;

    // Create dataspace for the field trace timing.
    hsize_t t0dims[2];
    t0dims[0] = Nx;
    t0dims[1] = Ny;
    dataspace = H5Screate_simple (2, t0dims, NULL);
    if (dataspace<0) throw(Screen_FileWriteError());
    // buffer the data
    buffer = new double[Nx*Ny];
    for (int ix=0; ix<Nx; ix++)
	    for (int iy=0; iy<Ny; iy++)
	        buffer[ix*Ny+iy] = A[ix][iy].get_t0();
    // Create the dataset creation property list
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(Screen_FileWriteError());
    // Create the dataset.
    dataset = H5Dcreate(file,
        "ObservationTime",		// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        dataspace, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dataset<0) throw(Screen_FileWriteError());
    // Write the data to the dataset
    status = H5Dwrite (dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        dataspace,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(Screen_FileWriteError());
    // attach scalar attributes
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "Nt", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    int Nt = get_Nt();
    status = H5Awrite(att, H5T_NATIVE_INT, &Nt);
    if (status<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "dt", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    double dt = get_dt();
    status = H5Awrite(att, H5T_NATIVE_DOUBLE, &dt);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Sclose (atts);
    if (status<0) throw(Screen_FileWriteError());
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Dclose (dataset);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Sclose (dataspace);
    if (status<0) throw(Screen_FileWriteError());
    delete buffer;

    // Create dataspace for the field data.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t dims[4];
    dims[0] = Nx;
    dims[1] = Ny;
    dims[2] = A[0][0].get_N();
    dims[3] = 6;
    dataspace = H5Screate_simple (4, dims, NULL);
    if (dataspace<0) throw(Screen_FileWriteError());
    // buffer the data
    buffer = new double[getBufferSize()];
    if (buffer==0) throw(Screen_MemoryAllocationError());
    bufferArray(buffer);
    // Create the dataset creation property list
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(Screen_FileWriteError());
    // Create the dataset.
    dataset = H5Dcreate (file,
        "ElMagField", 			// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        dataspace, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dataset<0) throw(Screen_FileWriteError());
    // Write the data to the dataset
    status = H5Dwrite (dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        dataspace,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(Screen_FileWriteError());
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Dclose (dataset);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Sclose (dataspace);
    if (status<0) throw(Screen_FileWriteError());
    delete buffer;

    status = H5Fclose (file);
    if (status<0) throw(Screen_FileWriteError());
    // no errors have occured if we made it 'til here
    cout << "writing HDF5 done." << endl;
    return;
}

