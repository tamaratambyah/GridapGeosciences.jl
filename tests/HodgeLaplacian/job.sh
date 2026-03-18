#!/bin/bash
#PBS -P bt62
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l ncpus=8
#PBS -l mem=32gb
#PBS -N hodge
#PBS -l wd

source $HOME/scripts/load-configs.sh
source $HOME/scripts/load-intel.sh

mpiexec -n $PBS_NCPUS julia --project=$PBS_O_WORKDIR -e'
  using MPI
  using PartitionedArrays
  include("debug.jl")
    
' 