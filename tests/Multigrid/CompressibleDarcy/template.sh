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

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

julia --project={{{projectdir}}} -O3 -e\


mpiexec -n {{ncpus}} julia --project={{{projectdir}}} -O3 -e\
  '
  using MPI

 include({{{driver}}})

  MPI.Init()
  np = MPI.Comm_size(MPI.COMM_WORLD)
  ranks = distribute_with_mpi(LinearIndices((np,)))
  
  MPI.Barrier(MPI.COMM_WORLD)

  convergence(ranks;
      c = {{c}},
      α = {{α}},
      n = {{n}},
      order = {{order}},
      iters = {{iters}},
      itu = {{itu}},
      itp = {{itp}},
      dir = {{{dir}}},
      return_vtk = {{return_vtk}},
      simName = "{{simName}}",
    )


  '
