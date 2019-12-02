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

DONE:
FieldTrace::*
void Screen::propagate_to(Vector target_pos, FieldTrace *target_trace)
void Screen::writeFieldHDF5(std::string filename)
GaussianWavePacket

TODO:
Screen::Screen(std::string filename)
FieldTrace Screen::dx_A(int ix, int iy)
FieldTrace Screen::dy_A(int ix, int iy)
void Screen::computeDerivatives()
TiltedScreen
PropagateStraight
plot_PowerDensity.py
plot_Screen_TD.py
TimeDomainField.py

