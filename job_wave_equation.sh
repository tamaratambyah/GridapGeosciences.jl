#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=8:00:00
#PBS -l ncpus=6
#PBS -l mem=24gb
#PBS -N wave_equation
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 6 julia --project=$PBS_O_WORKDIR -e'
  using MPI
  using PartitionedArrays
  include("tests/Geophysical/TransientWaveEquation.jl")
 
   with_mpi() do distribute        
        TransientWaveEquation.main_transient(distribute,6;n_ref_lvls=6,return_vtk=true)
  end                  

' 
