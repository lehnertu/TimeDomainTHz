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
    tri_ref empty;
    empty.p1=-1; empty.p2=-1; empty.p2=-1;
    triangles = std::vector<tri_ref>(Np,empty);
    field_points = std::vector<Vector>(Np);
    t0 = std::vector<double>(Nt,0.0);
    FieldTrace zero_trace(0.0,dt,Nt);
    A = std::vector<FieldTrace*>(Np);
    for (int i=0; i<Np; i++)
        A[i] = new FieldTrace(zero_trace);
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
    status = H5Aclose(attr);
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
    
    // read the position dataset
    dataset = H5Dopen2(file, "ObservationPosition", H5P_DEFAULT);
    if (dataset<0) throw(Screen_FileReadError());
    attr = H5Aopen_name(dataset, "Np");
    if (attr<0) throw(Screen_FileReadError());
    status =  H5Aread(attr, H5T_NATIVE_INT, &Np);
    if (status<0) throw(Screen_FileReadError());
    status = H5Aclose(attr);
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

    // read the mesh triangles
    dataset = H5Dopen2(file, "MeshTriangles", H5P_DEFAULT);
    if (dataset<0) throw(Screen_FileReadError());
    int *buf = new int[3*Np];
    status = H5Dread (dataset,
        H5T_NATIVE_INT, 		// mem type id
        H5S_ALL, 			    // mem space id
        H5S_ALL,
        H5P_DEFAULT,			// data transfer properties
        buf);
    status = H5Dclose(dataset);
    if (status<0) throw(Screen_FileReadError());
    triangles = std::vector<tri_ref>(Np);
    int *bt = buf;
    for (int ip=0; ip<Np; ip++)
    {   
        tri_ref tri;
        tri.p1 = *bt++;
        tri.p2 = *bt++;
        tri.p3 = *bt++;
        triangles[ip] = tri;
    }
    delete buf;
    
    // read the trace timing dataset
    // it is not an error if there are no data available
    if (H5Lexists(file, "ObservationTime", H5P_DEFAULT)<=0)
    {
        // no data available
        if (DEBUGLEVEL>=1) cout << "dataset ObservationTime not existing" << std::endl;
        Nt = 0;
        dt = 0.0;
        t0 = std::vector<double>(Np,0.0);
    }
    else
    {
        dataset = H5Dopen2(file, "ObservationTime", H5P_DEFAULT);
        attr = H5Aopen_name(dataset, "Nt");
        if (attr<0) throw(Screen_FileReadError());
        status =  H5Aread(attr, H5T_NATIVE_INT, &Nt);
        if (status<0) throw(Screen_FileReadError());
        status = H5Aclose(attr);
        if (status<0) throw(Screen_FileReadError());
        attr = H5Aopen_name(dataset, "dt");
        if (attr<0) throw(Screen_FileReadError());
        status =  H5Aread(attr, H5T_NATIVE_DOUBLE, &dt);
        if (status<0) throw(Screen_FileReadError());
        status = H5Aclose(attr);
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
    if (H5Lexists(file, "ElMagField", H5P_DEFAULT)<=0)
    {
        // no data available
        if (DEBUGLEVEL>=1) cout << "dataset ElMagField not existing" << std::endl;
        FieldTrace zero_trace(0.0,dt,Nt);
        A = std::vector<FieldTrace*>(Np);
        for (int i=0; i<Np; i++)
            A[i] = new FieldTrace(zero_trace);
    }
    else
    {
        dataset = H5Dopen2(file, "ElMagField", H5P_DEFAULT);
        ElMagField *field = new ElMagField[Np*Nt];
        status = H5Dread (dataset,
            H5T_NATIVE_DOUBLE, 		// mem type id
            H5S_ALL, 			    // mem space id
            H5S_ALL,
            H5P_DEFAULT,			// data transfer properties
            field);
        status = H5Dclose(dataset);
        if (status<0) throw(Screen_FileReadError());
        // transfer the data into the internal data structures
        FieldTrace trace(0.0, dt, Nt);
        A = std::vector<FieldTrace*>(Np);
        ElMagField *buf = field;
        for (int ip=0; ip<Np; ip++)
        {   
            trace.set_t0(t0[ip]);
            trace.set(buf,Nt);
            buf += Nt;
            A[ip] = new FieldTrace(trace);
        }
        delete field;
    }
        
    status = H5Fclose (file);
    if (status<0) throw(Screen_FileReadError());
}

void Screen::init()
{
    area = std::vector<double>(Np);
    total_area = 0.0;
    normal = std::vector<Vector>(Np);
    xi = std::vector<Vector>(Np);
    eta = std::vector<Vector>(Np);
    Vector avg_normal = Vector(0.0,0.0,0.0);
    for (int ip=0; ip<Np; ip++)
    {
        tri_ref tri = triangles[ip];
        Vector p1 = triangle_points[tri.p1];
        Vector p2 = triangle_points[tri.p2];
        Vector p3 = triangle_points[tri.p3];
        Vector n = cross( (p3-p1), (p2-p1));
        area[ip] = n.norm()*0.5;
        total_area += area[ip];
        n.normalize();
        normal[ip] = n;
        avg_normal += n;
    };
    avg_normal.normalize();
    // check all normals are angled less than 45 deg from the average
    for (int ip=0; ip<Np; ip++)
    {
        if (dot(normal[ip],avg_normal)<0.7) throw(Screen_NormalsAlignmentError());
    };
    // setup the cell-local coordinate systems
    if (avg_normal.y < 0.7)
    {
        // if the average normal vector is sufficiently far off the y direction
        // xi is confined to lay in the x-z-plane
        for (int ip=0; ip<Np; ip++)
        {
            xi[ip] = cross(Vector(0.0,1.0,0.0),normal[ip]);
            xi[ip].normalize();
            eta[ip] = cross(normal[ip],xi[ip]);
            eta[ip].normalize();
        }
    }
    else
    {
        // otherwise eta is confined to lay in the y-z-plane
        for (int ip=0; ip<Np; ip++)
        {
            eta[ip] = cross(normal[ip],Vector(1.0,0.0,0.0));
            eta[ip].normalize();
            xi[ip] = cross(eta[ip],normal[ip]);
            xi[ip].normalize();
        }
    }
}

Screen::~Screen()
{
    for (int i=0; i<Np; i++) delete A[i];
}

int Screen::get_Neighbourhood(int index, int *buffer)
{
    // the indexed triangle
    int i1 = triangles[index].p1;
    int i2 = triangles[index].p2;
    int i3 = triangles[index].p3;
    int count=0;
    int *buf = buffer;
    for (int i=0; i<Np; i++)
    {
        // test triangle
        int ti1 = triangles[i].p1;
        int ti2 = triangles[i].p2;
        int ti3 = triangles[i].p3;
        // if there are any common points
        if ( (ti1==i1) || (ti1==i2) || (ti1==i3) ||
             (ti2==i1) || (ti2==i2) || (ti2==i3) ||
             (ti3==i1) || (ti3==i2) || (ti3==i3) )
        {
            *buf++ = i;
            count++;
        }
    }
    return count;
}

double Screen::totalEnergy()
{
    double total = 0.0;
    for (int ip=0; ip<Np; ip++)
    {
        Vector S = A[ip]->Poynting();
        total -= area[ip]*dot(S,normal[ip]);
    }
    return total;
}

void Screen::writeReport(std::ostream *st)
{
    *st << "TimeDomainTHz - Screen" << std::endl;
    *st << "----------------------" << std::endl;
    *st << "Np=" << Np << "  Nt=" << get_Nt() << std::endl;
    st->precision(4);
    // *st << "average normal = (" << avg_normal.x << ", " << avg_normal.y << ", " << avg_normal.z << ")" << std::endl;
    st->precision(6);
    *st << "total area = " << total_area*1.0e4 << " cm²" << std::endl;
    *st << "Energy incident on screen = " << totalEnergy()*1.0e6 << " µJ" << std::endl;
    // find the cell with the highest power density
    int peak_index=-1;
    double peak=0.0;
    Vector Sp;
    for (int ip=0; ip<Np; ip++)
    {   
        Vector p = A[ip]->Poynting();
        if (p.norm()>peak)
        {
            Sp = p;
            peak = p.norm();
            peak_index = ip;
        }
    };
    st->precision(6);
    if (peak_index >= 0)
    {
        Vector p = get_point(peak_index);
        *st << "peak pos["  << peak_index << "] =   (" << p.x << ", " << p.y << ", " << p.z << ") m   " << std::endl;
        *st << "peak energy density = " << peak << " J/m²" << std::endl;
        *st << "Poynting vector S = (" << Sp.x << ", " << Sp.y << ", " << Sp.z << ") J/m²" << std::endl;
        *st << "t0  = " << A[peak_index]->get_t0()*1e9 << " ns  " << std::endl;
        *st << "dt  = " << dt*1.0e9 << " ns" << std::endl;
    };
    *st << std::endl;
}

void Screen::writeTraceReport(std::ostream *st, int ip)
{
    if ((ip>=0) and (ip<Np))
    {
        *st << "trace[" << ip << "]  ";
        Vector p = get_point(ip);
        *st << "pos=(" << p.x << ", " << p.y << ", " << p.z << ") m   " << std::endl;
        *st << "t0=" << A[ip]->get_t0()*1e9 << " ns  ";
        *st << "P=" << A[ip]->Poynting().norm() << " J/m²" << std::endl;
    }
}

void Screen::writeFieldHDF5(std::string filename)
{
    herr_t status;
    hid_t atts, att;
    hid_t dataspace, dcpl, dataset;
    
    cout << "writing HDF5 file " << filename << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // hid_t file = H5Fopen (filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file<0) throw(Screen_FileWriteError());

    // write the mesh corner points
    if (DEBUGLEVEL>=2) cout << "writing mesh points" << endl;
    hsize_t cdims[2];
    cdims[0] = Ncp;
    cdims[1] = 3;
    // create dataspace
    // Setting maximum size to NULL sets the maximum size to be the current size.
    dataspace = H5Screate_simple (2, cdims, NULL);
    if (dataspace<0) throw(Screen_FileWriteError());
    // buffer the data
    double *cbuffer = new double[Ncp*3];
    double *cbp = cbuffer;
    for (int icp=0; icp<Ncp; icp++)
	{
        Vector pos = triangle_points[icp];
        *cbp++ = pos.x;
        *cbp++ = pos.y;
        *cbp++ = pos.z;
    };
    // Create the dataset creation property list
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(Screen_FileWriteError());
    // Create the dataset.
    dataset = H5Dcreate(file,
        "MeshCornerPoints",		// dataset name
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
        cbuffer);
    if (status<0) throw(Screen_FileWriteError());
    // attach scalar attributes
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "Ncp", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    status = H5Awrite(att, H5T_NATIVE_INT, &Ncp);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Aclose (att);
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
    delete cbuffer;

    // write the position dataset
    if (DEBUGLEVEL>=2) cout << "writing mesh center points" << endl;
    hsize_t pdims[2];
    pdims[0] = Np;
    pdims[1] = 3;
    // create dataspace
    // Setting maximum size to NULL sets the maximum size to be the current size.
    dataspace = H5Screate_simple (2, pdims, NULL);
    if (dataspace<0) throw(Screen_FileWriteError());
    // buffer the data
    double *pbuffer = new double[Np*3];
    double *pbp = pbuffer;
    for (int ip=0; ip<Np; ip++)
	{
        Vector pos = field_points[ip];
        *pbp++ = pos.x;
        *pbp++ = pos.y;
        *pbp++ = pos.z;
    };
    // Create the dataset creation property list
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(Screen_FileWriteError());
    // Create the dataset.
    dataset = H5Dcreate(file,
        "ObservationPosition",	// dataset name
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
        pbuffer);
    if (status<0) throw(Screen_FileWriteError());
    // attach scalar attributes
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "Np", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    status = H5Awrite(att, H5T_NATIVE_INT, &Np);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Aclose (att);
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
    delete pbuffer;

    // write the mesh triangles
    if (DEBUGLEVEL>=2) cout << "writing mesh triangles" << endl;
    hsize_t tdims[2];
    tdims[0] = Np;
    tdims[1] = 3;
    // create dataspace
    // Setting maximum size to NULL sets the maximum size to be the current size.
    dataspace = H5Screate_simple (2, tdims, NULL);
    if (dataspace<0) throw(Screen_FileWriteError());
    // buffer the data
    int *tribuffer = new int[Np*3];
    int *tribp = tribuffer;
    for (int ip=0; ip<Np; ip++)
	{
	    tri_ref tri = triangles[ip];
        *tribp++ = tri.p1;
        *tribp++ = tri.p2;
        *tribp++ = tri.p3;
    };
    // Create the dataset creation property list
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(Screen_FileWriteError());
    // Create the dataset.
    dataset = H5Dcreate(file,
        "MeshTriangles",	    // dataset name
        H5T_NATIVE_INT, 		// data type
        dataspace, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dataset<0) throw(Screen_FileWriteError());
    // Write the data to the dataset
    status = H5Dwrite (dataset,
        H5T_NATIVE_INT, 		// mem type id
        H5S_ALL, 			    // mem space id
        dataspace,
        H5P_DEFAULT,			// data transfer properties
        tribuffer);
    if (status<0) throw(Screen_FileWriteError());
    // attach scalar attributes
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "Ntri", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    status = H5Awrite(att, H5T_NATIVE_INT, &Np);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Aclose (att);
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
    delete tribuffer;

    // write the timing dataset
    if (DEBUGLEVEL>=2) cout << "writing timing" << endl;
    hsize_t t0dims[1];
    t0dims[0] = Np;
    // create dataspace
    // Setting maximum size to NULL sets the maximum size to be the current size.
    dataspace = H5Screate_simple (1, t0dims, NULL);
    if (dataspace<0) throw(Screen_FileWriteError());
    // buffer the data
    double *t0buffer = new double[Np];
    double *t0bp = t0buffer;
    for (int ip=0; ip<Np; ip++)
        *t0bp++ = t0[ip];
    // Create the dataset creation property list
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(Screen_FileWriteError());
    // Create the dataset.
    dataset = H5Dcreate(file,
        "ObservationTime",	    // dataset name
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
        t0buffer);
    if (status<0) throw(Screen_FileWriteError());
    // attach scalar attributes
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "Nt", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    status = H5Awrite(att, H5T_NATIVE_INT, &Nt);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Aclose (att);
    if (status<0) throw(Screen_FileWriteError());
    att = H5Acreate2(dataset, "dt", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (att<0) throw(Screen_FileWriteError());
    status = H5Awrite(att, H5T_NATIVE_DOUBLE, &dt);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Aclose (att);
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
    delete t0buffer;

    // write the fields dataset
    if (DEBUGLEVEL>=2) cout << "writing fields" << endl;
    hsize_t fdims[3];
    fdims[0] = Np;
    fdims[1] = Nt;
    fdims[2] = 6;
    // create dataspace
    // Setting maximum size to NULL sets the maximum size to be the current size.
    dataspace = H5Screate_simple (3, fdims, NULL);
    if (dataspace<0) throw(Screen_FileWriteError());
    // buffer the data
    double *fbuffer = new double[Np*Nt*6];
    if (fbuffer==0) throw(Screen_MemoryAllocationError());
    double *fbp = fbuffer;
    for (int ip=0; ip<Np; ip++)
    {
        FieldTrace trace = get_trace(ip);
        for (int it=0; it<Nt; it++)
        {
            ElMagField f = trace.get_field(it);
            *fbp++ = f.E().x;
            *fbp++ = f.E().y;
            *fbp++ = f.E().z;
            *fbp++ = f.B().x;
            *fbp++ = f.B().y;
            *fbp++ = f.B().z;
        }
    }
    // Create the dataset creation property list
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(Screen_FileWriteError());
    // Create the dataset.
    dataset = H5Dcreate(file,
        "ElMagField",           // dataset name
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
        fbuffer);
    if (status<0) throw(Screen_FileWriteError());
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Dclose (dataset);
    if (status<0) throw(Screen_FileWriteError());
    status = H5Sclose (dataspace);
    if (status<0) throw(Screen_FileWriteError());
    delete fbuffer;

    // we are done with the file
    status = H5Fclose (file);
    if (status<0) throw(Screen_FileWriteError());

    // no errors have occured if we made it 'til here
    cout << "writing HDF5 done." << endl;
    return;
}

