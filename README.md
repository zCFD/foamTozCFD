# foamTozCFD
OpenFOAM to zCFD converter

## Build Dependencies:

hdf5 from https://support.hdfgroup.org/HDF5/ 

set environment variable HDF5_HOME to install location

hdf5 helper from https://github.com/zenotech/hdf5

Install/clone in $HDF5_HOME/include

## Build instructions:

source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc

wmake

## Run instructions:

source ~/OpenFOAM/OpenFOAM-X.X/etc/bashrc

foamTozCFD
