# foamTozCFD
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

### Mesh conversion

```
source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc
foamTozCFD

```

### Export solution


foamTozCFD
