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
- ghost pulse at the border of the screen (max 25 ps late) 

Analysis:

- first term is 4 orders of magnitude weaker than all 3 combined
- first term has about the same relative strength of the ghost pulse
- first term does not show the excessive longitudinal E component

- second term no excessive longitudinal E component
- third term excessive longitudinal E component !!!

- all terms show approximately equal relative strength of the ghost pulse
- ghost pulse is constant absolute strength over area
- gost pulse varies in delay depending on the position (latest far out)

=> the formula for the normal derivative is supposed to work on field
components in the local coordinate system - but global coordinates are used.

BUGFIX: compute spatial derivatives using screen-local coordinates for the field components

=>  excessive longitudinal components fixed
    ghost pulses still visible at the edges of the observation screens
    (likely diffraction artifacts, obscured in regions of greater intensity)
    intensities are correct
    BUG: polarity is reversed when comparing a single and a double propagation step

