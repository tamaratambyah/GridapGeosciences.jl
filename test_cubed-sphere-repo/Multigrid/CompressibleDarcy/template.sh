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

  include("{{{driver}}}")
  
  MPI.Barrier(MPI.COMM_WORLD)

  convergence(ranks;
      c = {{c}},
      α = {{α}},
      n = {{n}},
      order = {{order}},
      iters = {{iters}},
      itu = {{itu}},
      itp = {{itp}},
      dir = "{{{dir}}}",
      return_vtk = {{return_vtk}},
      simName = "{{simName}}",
    )

  i_am_main(ranks) && println("--DONE--")

'
