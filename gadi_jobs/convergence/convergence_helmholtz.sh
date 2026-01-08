#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l ncpus=4
#PBS -l mem=16gb
#PBS -N helmholtz
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 4 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/Laplace/Helmholtz.jl")

    with_mpi() do distribute
        Helmholtz.main(distribute,4;octree=true)
    end

' 
