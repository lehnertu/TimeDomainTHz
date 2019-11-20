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
}

Screen::~Screen()
{
}

Vector Screen::get_center()
{
    return Center;
}

Vector Screen::get_point(int ix, int iy)
{
    return Center + xVec*((double)ix-0.5*(double)Nx) + yVec*((double)iy-0.5*(double)Ny);
}

void Screen::set_Trace(int ix, int iy, FieldTrace trace)
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
            sum += dot(A[ix][iy].Poynting(),Normal);
    return(sum*dA);
}

void Screen::writeFieldHDF5()
{
    herr_t status;
    cout << "writing HDF5 file " << FileName << endl;
    // Create a new file using the default properties.
    hid_t file = H5Fcreate (FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Fcreate()"));
    
    // Create dataspace for the observation positions.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t pdims[3];
    pdims[0] = Nx;
    pdims[1] = Ny;
    pdims[2] = 3;
    hid_t pspace = H5Screate_simple (3, pdims, NULL);
    if (pspace<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(pspace)"));
    // buffer the data
    double *buffer = new double[Nx*Ny*3];
    double *bp = buffer;
    for (unsigned int ix = 0; ix < Nx; ix++)
	for (unsigned int iy = 0; iy < Ny; iy++)
	{
	    Vector pos = CellPosition(ix, iy);
	    *bp++ = pos.x;
	    *bp++ = pos.y;
	    *bp++ = pos.z;
	};
    // Create the dataset creation property list
    hid_t pdcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (pdcpl<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pcreate(pdcpl)"));
    // Create the dataset.
    hid_t pdset = H5Dcreate(file,
        "ObservationPosition",		// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        pspace, H5P_DEFAULT,
        pdcpl, H5P_DEFAULT);
    if (pdset<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dcreate(pdset)"));
    // Write the data to the dataset
    status = H5Dwrite (pdset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        pspace,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dwrite(pdset)"));
    // attach scalar attributes
    hid_t atts1  = H5Screate(H5S_SCALAR);
    if (atts1<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(atts1)"));
    hid_t attr1 = H5Acreate2(pdset, "Nx", H5T_NATIVE_INT, atts1, H5P_DEFAULT, H5P_DEFAULT);
    if (attr1<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(attr1)"));
    status = H5Awrite(attr1, H5T_NATIVE_INT, &Nx);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(attr1)"));
    status = H5Sclose (atts1);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(atts1)"));
    hid_t atts2  = H5Screate(H5S_SCALAR);
    if (atts2<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(atts2)"));
    hid_t attr2 = H5Acreate2(pdset, "Ny", H5T_NATIVE_INT, atts2, H5P_DEFAULT, H5P_DEFAULT);
    if (attr2<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(attr2)"));
    status = H5Awrite(attr2, H5T_NATIVE_INT, &Ny);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(attr2)"));
    status = H5Sclose (atts2);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(atts2)"));
    // Close and release resources.
    status = H5Pclose (pdcpl);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pclose(pdcpl)"));
    status = H5Dclose (pdset);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(pdset)"));
    status = H5Sclose (pspace);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(pspace)"));
    delete[] buffer;

    // Create dataspace for the field data.
    // Setting maximum size to NULL sets the maximum size to be the current size.
    hsize_t dims[4];
    dims[0] = Nx;
    dims[1] = Ny;
    dims[2] = NOTS;
    dims[3] = 6;
    hid_t space = H5Screate_simple (4, dims, NULL);
    if (space<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(space)"));
    // buffer the data
    buffer = getBuffer();
    // Create the dataset creation property list
    hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
    if (dcpl<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pcreate(dcpl)"));
    // Create the dataset.
    hid_t dset = H5Dcreate (file,
        "ElMagField", 			// dataset name
        H5T_NATIVE_DOUBLE,		// data type
        space, H5P_DEFAULT,
        dcpl, H5P_DEFAULT);
    if (dset<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dcreate(dset)"));
    // Write the data to the dataset
    status = H5Dwrite (dset,
        H5T_NATIVE_DOUBLE, 		// mem type id
        H5S_ALL, 			// mem space id
        space,
        H5P_DEFAULT,			// data transfer properties
        buffer);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dwrite(dset)"));
    // attach scalar attributes
    hid_t atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(t0)"));
    hid_t attr = H5Acreate2(dset, "t0", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(t0)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &t0_obs);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(t0)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(t0)"));
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(dt)"));
    attr = H5Acreate2(dset, "dt", H5T_NATIVE_DOUBLE, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(dt)"));
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dt_obs);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(dt)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(dt)"));
    atts  = H5Screate(H5S_SCALAR);
    if (atts<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Screate(NOTS)"));
    attr = H5Acreate2(dset, "NOTS", H5T_NATIVE_INT, atts, H5P_DEFAULT, H5P_DEFAULT);
    if (attr<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Acreate2(NOTS)"));
    status = H5Awrite(attr, H5T_NATIVE_INT, &NOTS);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Awrite(NOTS)"));
    status = H5Sclose (atts);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(NOTS)"));
    // Close and release resources.
    status = H5Pclose (dcpl);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Pclose(dcpl)"));
    status = H5Dclose (dset);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Dclose(dset)"));
    status = H5Sclose (space);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Sclose(space)"));
    delete[] buffer;

    status = H5Fclose (file);
    if (status<0) throw(IOexception("ScreenObserver::WriteTimeDomainFieldHDF5 - error in H5Fclose()"));
    // no errors have occured if we made it 'til here
    cout << "writing HDF5 done." << endl;
    return;
}

