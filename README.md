# transient-heatflow-moving-grid
Fortran code for transient heat-flow modelling with advective remapping, associated with Anees et al (2026), Solid Earth.
# Transient 1-D Lithospheric Heat-Flow Model

Fortran code for modelling transient one-dimensional conductive heat transport through a layered lithosphere.

The model solves the time-dependent heat equation using an implicit finite-difference formulation and a tridiagonal matrix solver. Thermal conductivity may vary with temperature, radiogenic heat production and material properties may vary between layers, and either fixed-temperature or basal heat-flow boundary conditions may be specified.

The code was developed by David Hindle, University of Göttingen.

## Repository contents

```text
.
├── README.md
├── Makefile
├── paramtrans.txt
├── transic.dat
├── src/
│   ├── pdtransient.f90
│   └── module_tdma.f90
└── plotting/
    └── pyplot.py
```

### `src/pdtransient.f90`

Main transient heat-flow model.

The program:

- reads the physical and numerical parameters from `paramtrans.txt`;
- reads the initial temperature profile from `transic.dat`;
- constructs the layered thermal-property model;
- solves the transient heat equation using an implicit finite-difference scheme;
- uses temperature-dependent thermal conductivity;
- permits either fixed-temperature or basal heat-flow lower boundary conditions;
- contains a coordinate-update formulation for vertical material motion;
- calculates temperature and heat-flow profiles through time;
- writes sequential model-output files.

### `src/module_tdma.f90`

Fortran module containing the tridiagonal matrix algorithm used to solve the implicit finite-difference system.

The implementation follows the tridiagonal methodology described by:

Sebben, S. and Baliga, B. R. (1995), *Some Extensions of Tridiagonal and Pentadiagonal Matrix Algorithms*, Numerical Heat Transfer, Part B, 28(3), 323–351.  
https://doi.org/10.1080/10407799508928837

### `paramtrans.txt`

Model parameter file.

The supplied example specifies:

- temperature-dependent/radiative thermal conductivity;
- grid spacing;
- timestep;
- model duration;
- number of saved solutions;
- surface and basal boundary conditions;
- vertical velocity;
- number of lithological layers;
- thermal conductivity;
- radiogenic heat production;
- density;
- specific heat capacity;
- depth of each layer boundary;
- properties of material entering the model domain.

Comments following the numerical values are ignored by the Fortran list-directed input.

### `transic.dat`

Initial temperature profile.

The supplied file contains one temperature value per numerical grid node. The example model contains 1501 nodes extending from the surface to 150 km depth at 100 m spacing.

Temperatures are given in kelvin.

### `plotting/pyplot.py`

Optional Python/Matplotlib script for plotting temperature and heat-flow profiles from selected output files.

This is not required to compile or run the Fortran model.

## Requirements

### Model

A Fortran compiler supporting Fortran 90 or later is required.

The code has been tested with:

```text
gfortran
```

GNU Make is also recommended.

### Plotting

The optional Python plotting script requires:

```text
Python 3
matplotlib
```

## Compilation

From the repository root, compile the model using:

```bash
make
```

This produces the executable:

```text
pdt.exe
```

Alternatively, the program can be compiled manually with `gfortran`.

For example:

```bash
gfortran -O3 -fbounds-check -o pdt.exe \
    src/module_tdma.f90 src/pdtransient.f90
```

The TDMA module must be compiled before the main program because `pdtransient.f90` contains:

```fortran
use TDMA
```

## Running the model

The executable expects the following files to be present in the working directory:

```text
paramtrans.txt
transic.dat
```

Run the model with:

```bash
./pdt.exe
```

The principal user-controlled parameters are contained in `paramtrans.txt`.

## Parameter file

The active part of `paramtrans.txt` has the following structure:

```text
radiative
bckind
kr
dx
dt
maxt
outno
bc1
bc2
v
number_of_layers

k  A  rho  Cp  depth
k  A  rho  Cp  depth
...
```

where:

- `radiative` switches the radiative contribution to temperature-dependent thermal conductivity on or off;
- `bckind = false` specifies a fixed-temperature lower boundary;
- `bckind = true` specifies a basal heat-flow boundary;
- `kr` is the temperature threshold for the radiative conductivity contribution, in K;
- `dx` is the nominal grid spacing, in metres;
- `dt` is the timestep, in years;
- `maxt` is the total model duration, in years;
- `outno` is the requested number of output intervals;
- `bc1` is the surface temperature boundary condition, in degrees Celsius;
- `bc2` is either basal temperature in degrees Celsius or basal heat flow, depending on `bckind`;
- `v` is vertical velocity in metres per year;
- `k` is thermal conductivity;
- `A` is volumetric radiogenic heat production;
- `rho` is density;
- `Cp` is specific heat capacity;
- `depth` is the depth to the base of the corresponding layer, in metres.

One additional property line following the defined layers specifies the properties assigned to material entering the model domain.

## Example model

The supplied parameter file defines four layers extending to 150 km depth with a grid spacing of 100 m.

The supplied initial-condition file therefore contains 1501 temperature values.

The supplied example has:

```text
dx       = 100 m
dt       = 100 yr
maxt     = 15 Myr
outno    = 10
Tsurface = 10 °C
Tbase    = 1400 °C
v        = 0 m yr^-1
```

Because `v = 0` in the supplied example, this particular configuration does not advect the grid.

## Moving-coordinate formulation

The program contains a time-dependent update of the vertical coordinate array `x`.

An older discrete-remapping section, in which temperature and material-property arrays are shifted between grid cells after sufficient accumulated displacement, is retained in the source code but is currently commented out.

The archived source should therefore be regarded as documenting the precise implementation contained in this repository rather than as a generic moving-mesh package.

## Output

At each requested save interval the program writes a file with names of the form:

```text
sol1001.dat
sol1002.dat
sol1003.dat
...
```

Each row contains four columns:

```text
depth_coordinate    temperature    thermal_conductivity    heat_flow
```

The program also writes:

```text
log.dat
```

containing diagnostic information about the grid and model setup.

Temperature is stored internally in kelvin.

## Plotting

The supplied Python script can be used to inspect selected output profiles:

```bash
python3 plotting/pyplot.py
```

Edit the `file_paths` list near the beginning of the script to select the output files to plot.

The first plot shows temperature versus depth and the second shows heat flow versus depth.

## Numerical method

The transient conductive heat equation is represented using an implicit finite-difference discretisation.

In schematic form,

```text
du/dt - d/dx(D du/dx) = q
```

is reduced at each timestep to a tridiagonal linear system.

The system is solved using the TDMA implementation in `module_tdma.f90`.

Thermal conductivity is recalculated as a function of temperature after each solution step.

## Boundary conditions

The surface boundary is a fixed-temperature (Dirichlet) condition.

The lower boundary may be either:

```text
bckind = false
```

for a fixed basal temperature, or

```text
bckind = true
```

for a specified basal heat-flow condition.

## Notes on reproducibility

The files `paramtrans.txt` and `transic.dat` supplied with this repository constitute the example input state archived with this version of the model.

For reproducibility, parameter files used for individual experiments should be retained together with the corresponding model version.

## Licence

See `LICENSE`.

## Citation

If you use this software, please cite the associated publication and the archived Zenodo release of this repository.

The permanent Zenodo DOI will be added here when the first repository release is archived.
