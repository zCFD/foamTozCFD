# foamTozCFD Converter
OpenFOAM to zCFD HDF5 converter

## Build Dependencies:

To build you will need a development build of OpenFoam and the following libraries:
- boost from https://boost.org/
- hdf5 from https://support.hdfgroup.org/HDF5/ 

Both of these are available from standard package repositories

Before compiling set environment variable HDF5_HOME to the install location of the hdf5 library.

## Build instructions:

```
source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc
wmake
```

## Run instructions:


To convert a case run foamTozCFD in your OpenFOAM case directory. The output will be created in a directory called "zCFDInterface" in the case directory. To get a list of the available options run.
```
source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc
foamTozCFD -help
```

### Mesh conversion


The mesh conversion will create a zCFD compatable HDF5 mesh file and also a python diction containing the boundary condition to zone mappings for zCFD.

```
source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc
foamTozCFD

```

### Export solution

The converter can also write the Foam solution out to zCFD format to allow for restarts from a Foam run. 

To convert the solution you will need to provide a reference velocity using the "-uinf" option. You can also supply a reference mach number (-minf) and density (-rinf), if you don't supply these then the values at sea level will be assumed. 

```
source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc
foamTozCFD -results -uinf 20

```

