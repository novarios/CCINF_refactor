#!/bin/bash
#
# build_andes.sh — Build CCDT on Andes (OLCF)
#
# Usage:
#   ./build_andes.sh              # clean + build (debug)
#   ./build_andes.sh release      # clean + build (release -O2)
#   ./build_andes.sh noclean      # build without cleaning (debug)
#   ./build_andes.sh clean        # clean only
#

module load gcc
module load openmpi
module load openblas
module load hdf5
module load gsl
export HDF5_DIR=${OLCF_HDF5_ROOT}

mkdir -p .mod

case "${1:-default}" in
    release)
        make -f makefile.andes clean
        make -f makefile.andes RELEASE=1
        ;;
    noclean)
        make -f makefile.andes
        ;;
    clean)
        make -f makefile.andes clean
        ;;
    *)
        make -f makefile.andes clean
        make -f makefile.andes
        ;;
esac
