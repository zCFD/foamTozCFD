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


To get a list of the available options run
```
source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc
foamTozCFD -help
```

### Mesh conversion

```
source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc
foamTozCFD

```

### Export solution

The converter can also write the Foam solution out to zCFD format to allow for restarts from a Foam run. 

```
source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc
foamTozCFD -results

```

