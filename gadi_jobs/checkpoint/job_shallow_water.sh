#!/bin/bash
#PBS -P bt62
#PBS -q normal
#PBS -l walltime=08:00:00
#PBS -l ncpus=96
#PBS -l mem=384gb
#PBS -N shallow_water_W5
#PBS -l wd

source $HOME/scripts/load-configs-bt62.sh
source $HOME/scripts/load-intel.sh

 
mpiexec -n $PBS_NCPUS julia --project=$PBS_O_WORKDIR -e'
  using PartitionedArrays
  using MPI

  MPI.Init()
  np = MPI.Comm_size(MPI.COMM_WORLD)
  ranks = distribute_with_mpi(LinearIndices((np,)))

  i_am_main(ranks) && println("--START--")

  include("test/TransientCheckingpointingTests/TransientShallowWater.jl") 

  with_mpi() do distribute
      main_transient(distribute,np;restart=false,
        n_ref_lvls=6,p_fe=1,return_vtk=1,freq=250) 
  end 

  i_am_main(ranks) && println("--DONE--")

'> /scratch/bt62/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/bt62/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
