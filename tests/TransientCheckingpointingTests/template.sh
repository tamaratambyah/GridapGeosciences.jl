#!/bin/bash
#PBS -P {{project}}
#PBS -q {{q}} 
#PBS -l walltime={{walltime}}
#PBS -l ncpus={{ncpus}}
#PBS -l mem={{mem}}
#PBS -N {{{name}}}
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh
 
mpiexec -n {{ncpus}} julia --project=$PBS_O_WORKDIR -e'
  using PartitionedArrays
  using MPI

  MPI.Init()
  np = MPI.Comm_size(MPI.COMM_WORLD)
  ranks = distribute_with_mpi(LinearIndices((np,)))

  i_am_main(ranks) && println("--START--")

  include("tests/TransientCheckingpointingTests/TransientThermalShallowWater.jl") 

  
  MPI.Barrier(MPI.COMM_WORLD)

    with_mpi() do distribute
        main_transient(distribute,{{ncpus}};restart=false,n_ref_lvls={{nl}},p_fe={{order}})
    end 

  i_am_main(ranks) && println("--DONE--")

' > {{{o}}}.out.${PBS_JOBID} 2> {{{e}}}.err.${PBS_JOBID}
