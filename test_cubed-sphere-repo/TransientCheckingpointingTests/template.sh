#!/bin/bash
#PBS -P {{project}}
#PBS -q {{q}} 
#PBS -l walltime={{walltime}}
#PBS -l ncpus={{ncpus}}
#PBS -l mem={{mem}}
#PBS -N {{{name}}}
#PBS -l wd
#PBS -o {{{o}}}
#PBS -e {{{e}}} 

source $HOME/scripts/load-configs-{{project}}.sh
source $HOME/scripts/load-intel.sh
 
mpiexec -n {{ncpus}} julia --project=$PBS_O_WORKDIR -e'
  using PartitionedArrays
  using MPI

  MPI.Init()
  np = MPI.Comm_size(MPI.COMM_WORLD)
  ranks = distribute_with_mpi(LinearIndices((np,)))

  i_am_main(ranks) && println("--START--")

  include("test/TransientCheckingpointingTests/TransientShallowWater.jl") 

  with_mpi() do distribute
      main_transient(distribute,np;restart=false,
        n_ref_lvls={{nl}},p_fe={{order}},return_vtk={{return_vtk}}) 
  end 

  i_am_main(ranks) && println("--DONE--")

' 
#> {{{o}}}.${PBS_JOBID} 2> {{{e}}}.${PBS_JOBID}
