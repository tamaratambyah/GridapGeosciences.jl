#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=24
#PBS -l mem=96gb
#PBS -N shallow_water
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 24 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/Geophysical/ShallowWater.jl")
    include("tests/Geophysical/TransientShallowWater.jl")

    with_mpi() do distribute
        # ShallowWater.main(distribute,4;octree=true)
        TransientShallowWater.main(distribute,24;octree=true)
    end

' 
