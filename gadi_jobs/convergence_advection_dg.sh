#!/bin/bash
#PBS -P bt62
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=24
#PBS -l mem=96gb
#PBS -N advection_dg
#PBS -l wd

source $HOME/scripts/load-configs.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 24 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/Advection/AdvectionDGUpwinding.jl")
    include("tests/Advection/TransientAdvectionDGUpwinding.jl")

    with_mpi() do distribute
        # AdvectionDGUpwinding.main(distribute,24;octree=true)
        TransientAdvectionDGUpwinding.main(distribute,24;octree=true)
    end

' 
