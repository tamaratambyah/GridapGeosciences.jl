#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=48
#PBS -l mem=192gb
#PBS -N shallow_water_latlon
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n $PBS_NCPUS julia --project=$PBS_O_WORKDIR -e'
  using MPI
  using PartitionedArrays
  include("tests/Geophysical/TransientShallowWater.jl")
 
   with_mpi() do distribute        
        TransientShallowWater.main_transient(distribute,48;octree=true,n_ref_lvls=6,return_vtk=true)
  end                  

' > /scratch/zg98/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/zg98/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
