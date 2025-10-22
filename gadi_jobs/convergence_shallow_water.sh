#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l ncpus=4
#PBS -l mem=16gb
#PBS -N shallow_water
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 4 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/Geophysical/ShallowWater.jl")
    include("tests/Geophysical/TransientShallowWater.jl")

    with_mpi() do distribute
        ShallowWater.main(distribute,4;octree=true)
        TransientShallowWater.main(distribute,4;octree=true)
    end

' 
