#!/bin/bash

module purge
module load Core/26.01
module load gcc/14.2.0
module load mpich/5.0.1
module load hdf5/1.14.6-mpi
module load openblas/0.3.33-omp
module load gsl/2.8

export HDF5_DIR=${OLCF_HDF5_ROOT}

noclean=false
if (( $# >= 1 )); then
        echo $1
        if [ "$1" == "noclean" ]; then
                noclean=true
                echo "noclean set -- will not run <<make clean>> first"
        fi
fi

if [ "$noclean" = false ]; then
        make -f makefile.riker clean
fi
make -f makefile.riker
