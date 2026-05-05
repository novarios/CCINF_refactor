#!/bin/bash
#
# build.sh — Build script for CCDT v2 code
#
# Usage:
#   ./build.sh              # build (debug mode)
#   ./build.sh release      # build with -O3 optimization
#   ./build.sh clean        # clean build artifacts
#

set -e

HDF5_DIR="/usr/lib/x86_64-linux-gnu/hdf5/openmpi"
HDF5_INCLUDEDIR="${HDF5_DIR}/include"
HDF5_LIBDIR="${HDF5_DIR}/lib"

NN_GRID_DIR="nn_grid_from_pw"

COMP="mpif90 -fopenmp -fdiagnostics-color=always -ffree-line-length-512 -fallow-argument-mismatch"
FCD="${COMP} -Og -g -fbounds-check -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
FC="${COMP} -O3"
XLF="${FCD}"

HDF5_INCLUDE="-I${HDF5_INCLUDEDIR}"
HDF5_LIBS="-L${HDF5_LIBDIR} ${HDF5_LIBDIR}/libhdf5hl_fortran.a ${HDF5_LIBDIR}/libhdf5_hl.a ${HDF5_LIBDIR}/libhdf5_fortran.a ${HDF5_LIBDIR}/libhdf5.a -ldl -lm -lz -Wl,-rpath -Wl,${HDF5_LIBDIR}"
LIBZ="/usr/lib/x86_64-linux-gnu/libz.a"
LIBSZ="/usr/lib/x86_64-linux-gnu/libsz.a"
LIBSA="/usr/lib/x86_64-linux-gnu/libaec.a"
LIB="${LIBZ} ${LIBSZ} ${LIBSA} -lrt -lm"
LOCAL_LIBS="-lopenblas -lgsl -lgslcblas"
LIBS="${LOCAL_LIBS} ${HDF5_LIBS} ${LIB}"

SOURCES=(
    modules.f90
    build_status.f90
    mem_tracker.f90
    contracts.f90
    mpi_check.f90
    hdf5_routines.f90
    chiral_module_andreas_sam.f90
    mappings.f90
    library.f90
    test.f90
    diis.f90
    paeom_routines.f90
    paeom_diagrams.f90
    preom_routines.f90
    preom_diagrams.f90
    eom_states.f90
    hbar_diagrams.f90
    ground_state_routines.f90
    mbpt.f90
    t3_diagrams.f90
    ground_state.f90
    energy.f90
    interaction.f90
    basis.f90
    ccdt_main.f90
)

do_clean() {
    echo "Cleaning..."
    rm -f *.o *.mod *.exe *__genmod.f90
    echo "Done."
}

do_build() {
    local mode="$1"
    if [ "$mode" = "release" ]; then
        XLF="${FC}"
        echo "=== CCDT v2 Build (RELEASE) ==="
    else
        echo "=== CCDT v2 Build (DEBUG) ==="
    fi
    echo "Compiler: $(mpif90 --version 2>&1 | head -1)"
    echo ""

    NN_OBJECTS=""
    if [ -d "${NN_GRID_DIR}" ]; then
        for f in ${NN_GRID_DIR}/*.f90 ${NN_GRID_DIR}/*.f; do
            [ -f "$f" ] || continue
            echo "  Compiling $(basename $f) ..."
            ${XLF} ${HDF5_INCLUDE} -I${NN_GRID_DIR} -c "$f"
        done
        NN_OBJECTS=$(ls ${NN_GRID_DIR}/*.o 2>/dev/null || echo "")
    fi

    for src in "${SOURCES[@]}"; do
        if [ ! -f "$src" ]; then
            echo "ERROR: Missing source file: $src"
            exit 1
        fi
        echo "  Compiling $src ..."
        ${XLF} ${HDF5_INCLUDE} -c "$src"
    done

    OBJECTS=""
    for src in "${SOURCES[@]}"; do
        OBJECTS="${OBJECTS} ${src%.f90}.o"
    done

    echo "  Linking ccdt.exe ..."
    ${XLF} -o ccdt.exe ${OBJECTS} ${NN_OBJECTS} ${LIBS}
    echo ""
    echo "=== Build successful: ./ccdt.exe ==="
}

case "${1:-build}" in
    release)  do_build release ;;
    clean)    do_clean ;;
    build|"") do_build debug ;;
    *)        echo "Usage: $0 [build|release|clean]"; exit 1 ;;
esac
