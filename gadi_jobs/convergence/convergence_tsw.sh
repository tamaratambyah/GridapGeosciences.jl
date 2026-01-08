#!/bin/bash
#PBS -P bt62
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l ncpus=24
#PBS -l mem=512gb
#PBS -N tsw
#PBS -l wd

source $HOME/scripts/load-configs.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 24 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays 
    include("tests/Geophysical/TransientThermalShallowWater.jl")

    with_mpi() do distribute
       main(distribute,24;octree=true)
    end

' 
