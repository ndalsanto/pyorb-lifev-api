#!/bin/bash

rm CMakeCache.txt
rm -rf CMakeFiles

# if you have installed LifeV by following the instructions on www.lifev.org, LIBRARIES_BASE_DIRECTORY is a sub-directory of your lifev-env folder

LIBRARIES_BASE_DIRECTORY=<base_directory_of_libraries>

LAPACK_LIBRARIES=${LIBRARIES_BASE_DIRECTORY}OpenBLAS-0.2.20_installRelease/lib/
METIS_INCLUDE_DIR=${LIBRARIES_BASE_DIRECTORY}metis-5.1.0_installRelease/include/
METIS_LIB_DIR=${LIBRARIES_BASE_DIRECTORY}metis-5.1.0_installRelease/lib/
PARMETIS_LIBRARIES=${LIBRARIES_BASE_DIRECTORY}parmetis-4.0.3_installRelease/lib/
PARMETIS_INCLUDE_DIR=${LIBRARIES_BASE_DIRECTORY}parmetis-4.0.3_installRelease/include/
PARMETIS_LIB_DIR=${LIBRARIES_BASE_DIRECTORY}parmetis-4.0.3_installRelease/lib/
HDF5_INCLUDE_DIR=${LIBRARIES_BASE_DIRECTORY}hdf5-1.8.19_installRelease/include/
HDF5_LIB_DIR=${LIBRARIES_BASE_DIRECTORY}hdf5-1.8.19_installRelease/lib/
TRILINOS_INCLUDE_DIR=${LIBRARIES_BASE_DIRECTORY}trilinos-release-12-12-1_installRelease/include/
TRILINOS_LIB_DIR=${LIBRARIES_BASE_DIRECTORY}trilinos-release-12-12-1_installRelease/lib/
MPI_INCLUDE_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}/openmpi/include/
MPI_LIB_DIRECTORY=${LIBRARIES_BASE_DIRECTORY}openmpi/lib/


LIFEV_INSTALLATION=<yout lifeV installation>
LIFEV_4_PYORB_DIR=<target build installation>

cmake \
-D TPL_LifeV_INCLUDE_DIRS:PATH=${LIFEV_INSTALLATION}include/ \
-D TPL_LifeV_LIBRARY_DIRS:PATH=${LIFEV_INSTALLATION}lib/ \
-D TPL_LAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES} \
-D TPL_METIS_INCLUDE_DIRS:PATH=${METIS_INCLUDE_DIR} \
-D TPL_METIS_LIBRARY_DIRS:PATH=${METIS_LIB_DIR} \
-D TPL_ParMETIS_LIBRARIES:PATH="${PARMETIS_LIBRARIES}" \
-D TPL_ParMETIS_INCLUDE_DIRS:PATH=${PARMETIS_INCLUDE_DIR} \
-D TPL_ParMETIS_LIBRARY_DIRS:PATH=${PARMETIS_LIB_DIR} \
-D TPL_HDF5_INCLUDE_DIRS:PATH=${HDF5_INCLUDE_DIR} \
-D TPL_HDF5_LIBRARY_DIRS:PATH=${HDF5_LIB_DIR} \
-D TPL_Trilinos_INCLUDE_DIRS:PATH=${TRILINOS_INCLUDE_DIR} \
-D TPL_Trilinos_LIBRARY_DIRS:PATH=${TRILINOS_LIB_DIR} \
-D TPL_MPI_INCLUDE_DIRS:PATH=${MPI_INCLUDE_DIRECTORY} \
-D TPL_MPI_LIBRARY_DIRS:PATH=${MPI_LIB_DIRECTORY} \
${LIFEV_4_PYORB_DIR}
