#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=36
#PBS -l mem=144gb
#PBS -N shallow_water
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 36 julia --project=$PBS_O_WORKDIR -e'
  using MPI
  using PartitionedArrays
  include("tests/Geophysical/TransientShallowWater.jl")
 
   with_mpi() do distribute        
        TransientShallowWater.main_transient(distribute,36;octree=true,n_ref_lvls=6,return_vtk=true)
  end                  

' 
