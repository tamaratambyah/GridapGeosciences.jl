#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l ncpus=6
#PBS -l mem=24gb
#PBS -N advection
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 6 julia --project=$PBS_O_WORKDIR -e'
  using MPI
  using PartitionedArrays
  options = """
    -ksp_type gmres
    -ksp_rtol 1.0e-6
    -ksp_converged_reason
    -ksp_monitor
    """

   with_debug() do distribute
      AdvectionSUPG.main(distribute;nprocs=6,options=options,
        n_ref_lvls=6,p_fe=1,CFL=0.1,tF=0.1,return_vtk=true)
  end                  

' 
