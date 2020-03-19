# TimeDomainTHz
THz radiation propagation in time-domain

The fields are defined on an unstructured triangular mesh.
The geometry is described by a list of corner points of the triangles
and a connection list indexing the 3 corner points of each triangle.
For every triangle a field trace is given as seen at the center point.

THz fields are described as time traces of (E,B) vector fields.
Every observation point has individual starting time for the trace.
This enables a great reduction in neccesary storage for screens
angled with respect to the direction of propagation.
All traces share a common time-step and length.

For propagation of radiation from one such mesh to another one
the Kirchhoff integral is evaluated directly.

file format
file.hdf5
- dataset "MeshCornerPoints"
    - attributes Ncp
    - array Ncp*3 double
- dataset "MeshTriangles"
    - attributes Ntri
    - array Ntri*3 uint32
- dataset "ObservationPosition"
    - attributes Np
    - array Np*3 double
- dataset "ObservationTime"
    - attributes Nt, dt
    - array Np double -- t0
- dataset "ElMagField"
    - no attributes 
    - array Np* Nt* 6 double -- FieldTrace[Np]

A few libraries are required to build the executables.

- For linear algebra calculations we use the [Eigen3](http://eigen.tuxfamily.org) library.
  It should preferably be installed through the systems package management system.
  In case this is not possible (i.e. no root access) it can be used from the git repository.
  ```cd lib```
  ```git clone https://gitlab.com/libeigen/eigen.git```
  A cmake script is provided which will find the library in either case.

- On most systems the [HDF5 library](https://support.hdfgroup.org/HDF5/) can be installed
  from the distribution repositories.
  It is necessary to set an environment variable for the library to be found by the build system

  ```export HDF5_ROOT=_path_to_library```


