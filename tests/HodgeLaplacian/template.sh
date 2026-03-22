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
 
mpiexec -n 1 julia --project=$PBS_O_WORKDIR -e'
  using PartitionedArrays
  using MPI

  MPI.Init()
  np = MPI.Comm_size(MPI.COMM_WORLD)
  ranks = distribute_with_mpi(LinearIndices((np,)))

  i_am_main(ranks) && println("--START--")

  include("{{{driver}}}")
  
  MPI.Barrier(MPI.COMM_WORLD)

  launch_hodge_laplacian(ranks,{{n}},
  {{order}},"{{{dir}}}",{{return_vtk}})

 
  i_am_main(ranks) && println("--DONE--")

'
