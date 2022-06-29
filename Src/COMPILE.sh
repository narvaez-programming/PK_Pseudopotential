#!/bin/bash

CC="gfortran"

FUNCTIONS="factorial_double.o factorial.o laplacian.o gaussian_normalization.o neumaier.o grid_rectangular.o grid_cubic.o convolve.o"
READ="read_file.o read_basis.o read_coords.o read_overlap_matrix.o read_kinetic_matrix.o read_mos.o read_pk_guess.o"
BUILD="construct_projection_matrix.o"
GRID_EVAL="core_mos_print.o pk.o pk_convolved.o pk_convolved_tapered.o"
TAPER='prep_tapering.o steinhauser.o taper_atom.o'

LAPACK="-llapack ${LAPACK_LIB}"
BLAS="-lblas ${BLAS_LIB}"

COMPILE="$CC -o pk.exe write_file.o constants.o $FUNCTIONS setup_variables.o $READ $BUILD find_pk_coefficients.o $GRID_EVAL $TAPER main.o $LAPACK $BLAS"

$CC -c write_file.f08 constants.f08
$CC -c grid_rectangular.f08 grid_cubic.f08
$CC -c factorial.f08 factorial_double.f08 neumaier.f08 laplacian.f08 gaussian_normalization.f08 convolve.f08
$CC -c setup_variables.f08
$CC -c read_file.f08 read_basis.f08 read_coords.f08 read_overlap_matrix.f08 read_kinetic_matrix.f08 read_mos.f08 read_pk_guess.f08
$CC -c construct_projection_matrix.f08
$CC -c find_pk_coefficients.f08
$CC -c core_mos_print.f08
$CC -c prep_tapering.f08 steinhauser.f08 taper_atom.f08
$CC -c pk.f08 pk_convolved.f08 pk_convolved_tapered.f08
$CC -c main.f08
  
$COMPILE
rm *.o *.mod 
