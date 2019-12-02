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
    // now we know the size of the dataset - create and fill a buffer
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
    Normal.normalize();
    Center = (pos11+pos00)*0.5;
    delete pos;

    // open the field dataset
    hid_t field_dataset = H5Dopen2(file, "ElMagField", H5P_DEFAULT);
    if (field_dataset<0) throw(Screen_FileReadError());
    int Nt;
    attr = H5Aopen_name(field_dataset, "NOTS");
    if (attr<0) throw(Screen_FileReadError());
    status =  H5Aread(attr, H5T_NATIVE_INT, &Nt);
    if (status<0) throw(Screen_FileReadError());
    double t0;
    attr = H5Aopen_name(field_dataset, "t0");
    if (attr<0) throw(Screen_FileReadError());
    status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &t0);
    if (status<0) throw(Screen_FileReadError());
    double dt;
    attr = H5Aopen_name(field_dataset, "dt");
    if (attr<0) throw(Screen_FileReadError());
    status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &dt);
    if (status<0) throw(Screen_FileReadError());
    // now we know the size of the dataset - create and fill a buffer
    ElMagField *field = new ElMagField[Nx*Ny*Nt];
    status = H5Dread (field_dataset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        field);
    // transfer the data into the internal data structures
    FieldTrace trace(t0, dt, Nt);
    std::vector<FieldTrace> zero_column(Ny, trace);
    A = std::vector< std::vector<FieldTrace> >(Nx,zero_column);
    ElMagField *buf = field;
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {   
            trace.set(buf,Nt);
            buf += Nt;
            A[ix][iy] = trace;
        }
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
    *st << "e_x = (" << xVec.x*1.0e3 << ", " << xVec.y*1.0e3 << ", " << xVec.z*1.0e3 << ") mm" << std::endl;
    *st << "e_y = (" << yVec.x*1.0e3 << ", " << yVec.y*1.0e3 << ", " << yVec.z*1.0e3 << ") mm" << std::endl;
    *st << "n   = (" << Normal.x << ", " << Normal.y << ", " << Normal.z << ")" << std::endl;
    st->precision(6);
    *st << "dA  = " << dA*1.0e6 << " mm²" << std::endl;
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
    FieldTrace trace = A[ix][iy];
    if (ix==0)
        trace = (A[1][iy] - A[0][iy]) / dX;
    else 
    {
        if (ix==Nx-1)
            trace = (A[Nx-1][iy] - A[Nx-2][iy]) / dX;
        else
            trace = (A[ix+1][iy] - A[ix-1][iy]) / (dX*2.0);
    }
    return trace;
}

FieldTrace Screen::dy_A(int ix, int iy)
{
    if (ix<0 || ix>=Nx) throw(Screen_IndexOutOfRange());
    if (iy<0 || iy>=Ny) throw(Screen_IndexOutOfRange());
    double dY = yVec.norm();
    FieldTrace trace = A[ix][iy];
    if (iy==0)
        trace = (A[ix][1] - A[ix][0]) / dY;
    else 
    {
        if (iy==Ny-1)
            trace = (A[ix][Ny-1] - A[ix][Ny-2]) / dY;
        else
            trace = (A[ix][iy+1] - A[ix][iy-1]) / (dY*2.0);
    }
    return trace;
}

void Screen::computeDerivatives()
{
    // compute the time-derivatives for all traces
    // pragma omp parallel for
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
            dt_A[ix][iy] = A[ix][iy].derivative();
    // compute the normal derivatives using the Maxwell equations
    // this is the drivative along the Normal direction (in the screen-local frame)
    // is is assumed that (xVec, yVec, Normal) form a right-handed orthogonal reference frame
    int Nt = A[0][0].get_N();
    FieldTrace trace(0.0,0.0,Nt);
    double *pEx, *pEy, *pEz, *pBx, *pBy, *pBz;
    double *dx_Ex = new double[Nx*Ny*Nt];
    pEx = dx_Ex;
    double *dx_Ez = new double[Nx*Ny*Nt];
    pEz = dx_Ez;
    double *dx_Bx = new double[Nx*Ny*Nt];
    pBx = dx_Bx;
    double *dx_Bz = new double[Nx*Ny*Nt];
    pBz = dx_Bz;
    // this is not parallelizable because the pointers are used by all threads
    // pragma omp parallel for private(trace)
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            trace = dx_A(ix,iy);
            for (int it=0; it<Nt; it++)
            {
                ElMagField f = trace.get_field(it);
                *pEx++ = f.E().x;
                *pEz++ = f.E().z;
                *pBx++ = f.B().x;
                *pBz++ = f.B().z;
            }
        };
    double *dy_Ey = new double[Nx*Ny*Nt];
    pEy = dy_Ey;
    double *dy_Ez = new double[Nx*Ny*Nt];
    pEz = dy_Ez;
    double *dy_By = new double[Nx*Ny*Nt];
    pBy = dy_By;
    double *dy_Bz = new double[Nx*Ny*Nt];
    pBz = dy_Bz;
    // this is not parallelizable because the pointers are used by all threads
    // pragma omp parallel for private(trace)
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            trace = dy_A(ix,iy);
            for (int it=0; it<Nt; it++)
            {
                ElMagField f = trace.get_field(it);
                *pEy++ = f.E().y;
                *pEz++ = f.E().z;
                *pBy++ = f.B().y;
                *pBz++ = f.B().z;
            }
        };
    double *dt_Ex = new double[Nx*Ny*Nt];
    pEx = dt_Ex;
    double *dt_Ey = new double[Nx*Ny*Nt];
    pEy = dt_Ey;
    double *dt_Bx = new double[Nx*Ny*Nt];
    pBx = dt_Bx;
    double *dt_By = new double[Nx*Ny*Nt];
    pBy = dt_By;
    // this is not parallelizable because the pointers are used by all threads
    // pragma omp parallel for private(trace)
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            trace = dt_A[ix][iy];
            for (int it=0; it<Nt; it++)
            {
                ElMagField f = trace.get_field(it);
                *pEx++ = f.E().x;
                *pEy++ = f.E().y;
                *pBx++ = f.B().x;
                *pBy++ = f.B().y;
            }
        };
    double *dz_Ex = new double[Nx*Ny*Nt];
    double *dz_Ey = new double[Nx*Ny*Nt];
    double *dz_Ez = new double[Nx*Ny*Nt];
    double *dz_Bx = new double[Nx*Ny*Nt];
    double *dz_By = new double[Nx*Ny*Nt];
    double *dz_Bz = new double[Nx*Ny*Nt];
    double cSquared = SpeedOfLight*SpeedOfLight;
    // this is not parallelizable because the pointers are used by all threads
    // pragma omp parallel for private(trace)
    for (int i=0; i<Nx*Ny*Nt; i++)
    {
        dz_Ex[i] = dx_Ez[i] - dt_By[i];
        dz_Ey[i] = dy_Ez[i] + dt_Bx[i];
        dz_Ez[i] = -dx_Ex[i] - dy_Ey[i];
        dz_Bx[i] = dx_Bz[i] + dt_Ey[i]/cSquared;
        dz_By[i] = dy_Bz[i] - dt_Ex[i]/cSquared;
        dz_Bz[i] = -dx_Bx[i] - dy_By[i];
    };
    pEx = dz_Ex;
    pEy = dz_Ey;
    pEz = dz_Ez;
    pBx = dz_Bx;
    pBy = dz_By;
    pBz = dz_Bz;
    // this is not parallelizable because the pointers are used by all threads
    // pragma omp parallel for
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            for (int it=0; it<Nt; it++)
            {
                Vector E = Vector(*pEx++,*pEy++,*pEz++);
                Vector B = Vector(*pBx++,*pBy++,*pBz++);
                ElMagField f = ElMagField(E,B);
                trace.set(it,f);
            }
            dn_A[ix][iy] = trace*(-1.0);
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
// FieldTrace Screen::propagation(Vector target, double t0_p, int N_p)
{
    // this is the result sum
    FieldTrace trace(target_trace);
    trace.zero();
    // this is the component to be added
    FieldTrace t2(target_trace);
    for (int ix=0; ix<Nx; ix++)
        for (int iy=0; iy<Ny; iy++)
        {
            Vector source = get_point(ix,iy);
            Vector RVec = target_pos - source;
            double R = RVec.norm();
            double R2 = R*R;
            double R3 = R2*R;
            FieldTrace t1 = A[ix][iy];
            t1.retarded(R/SpeedOfLight, &t2);
            trace += t2 * (dot(RVec,Normal)/R3);
            t1 = dt_A[ix][iy];
            t1.retarded(R/SpeedOfLight, &t2);
            trace += t2 * (dot(RVec,Normal)/(R2*SpeedOfLight));
            t1 = dn_A[ix][iy];
            t1.retarded(R/SpeedOfLight, &t2);
            trace += t2 * (-1.0/R);
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

