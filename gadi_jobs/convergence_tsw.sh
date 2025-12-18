#!/bin/bash
#PBS -P bt62
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=1
#PBS -l mem=96gb
#PBS -N tsw
#PBS -l wd

source $HOME/scripts/load-configs.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 1 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays 
    include("tests/Geophysical/TransientThermalShallowWater.jl")

    with_mpi() do distribute
       main(distribute,1;octree=true)
    end

' 
