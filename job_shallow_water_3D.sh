#!/bin/bash
#PBS -P bt62
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=48
#PBS -l mem=192gb
#PBS -N shallow_water_3D
#PBS -l wd

source $HOME/scripts/load-configs-$PROJECT.sh
source $HOME/scripts/load-intel.sh

mpiexec -n $PBS_NCPUS julia --project=$PBS_O_WORKDIR -e'
  using MPI
  using PartitionedArrays
  include("tests/Geophysical/TransientShallowWater_3D.jl")
 
  with_mpi() do distribute        
    main_transient(distribute,48;n_ref_lvls=6)
  end                  

' > /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
