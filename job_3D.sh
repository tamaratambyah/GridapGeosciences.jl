#!/bin/bash
#PBS -P bt62
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l ncpus=96
#PBS -l mem=512gb
#PBS -N linear_boussineq_amg
#PBS -l wd

source $HOME/scripts/load-configs.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 96 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/Geophysical/TransientLinearBoussinesq.jl")
    
    with_mpi() do distribute
       main_transient(distribute,96;n_ref_lvls=4,p_fe=1,CFL=0.2,return_vtk=true)
    end      
' > /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
