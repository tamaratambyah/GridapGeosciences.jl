#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l ncpus=4
#PBS -l mem=16gb
#PBS -N advection_dg
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 4 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/Advection/AdvectionDGUpwinding.jl")
    include("tests/Advection/TransientAdvectionDGUpwinding.jl")

    with_mpi() do distribute
        AdvectionDGUpwinding.main(distribute,4;octree=true)
        TransientAdvectionDGUpwinding.main(distribute,4;octree=true)
    end

' 
