export PROJECT="bt62"
export PETSC_VERSION="3.22.4"
export P4EST_VERSION="2.8.5"

SCRATCH="/scratch/$PROJECT/$USER"
export JULIA_DEPOT_PATH="$SCRATCH/.julia"

module purge
module load pbs

module load intel-compiler-llvm/2025.0.4
module load intel-mpi/2021.14.1
module load intel-mkl/2025.0.1

# For some reason are not set to mpi
export CONFIGURE_FLAGS='--with-cc=mpiicc --with-cxx=mpiicpc --with-fc=mpiifx --with-f77=mpiifort --with-f90=mpiifort'
export CC=mpiicc
export CXX=mpiicpc
export CP=mpiicpc
export F77=mpiifort
export F90=mpiifort
export FC=mpiifx

export MPI_VERSION="intel-$INTEL_MPI_VERSION"
export JULIA_MPI_PATH=$INTEL_MPI_ROOT

export JULIA_PETSC_LIBRARY="$HOME/bin/petsc/$PETSC_VERSION-$MPI_VERSION/lib/libpetsc"
export P4EST_ROOT_DIR="$HOME/bin/p4est/$P4EST_VERSION-$MPI_VERSION"

module load julia/1.11.0
