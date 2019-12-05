# TimeDomainTHz
THz radiation propagation in time-domain

THz fields are described as time traces of (E,B) vector fields
arranged in a (rectangular) grid. For propagation the Kirchhoff integral
is evaluated directly.

This feature branch implements individual timing for then traces
comprising one dataset (screen). This enables a great reduction
in neccesary storage for screens angled with respect to
the direction of propagation. Still, all traces have common
time-step and length.

file format
file.hdf5
    - dataset "ObservationPosition"
        - array Nx*Ny*3 double -- Vector[Nx,Ny]
        - attributes Nx, Ny
    - dataset "ObservationTime"
        - array Nx*Ny double -- t0
        - attributes Nt, dt
    - dataset "ElMagField"
        - array Nx*Ny*Nt*6 double -- FieldTrace[Nx,Ny]
        - no attributes 

BUG:

- intensity too large (factor sqrt(2)) when propagating from angled screen
- large longitudinal E field when propagating from angled screen
