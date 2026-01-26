#!/bin/bash
#PBS -P zg98
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l ncpus=208
#PBS -l mem=1024gb
#PBS -N linear_boussineq_amg
#PBS -l wd

source $HOME/scripts/load-configs.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 204 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays 
    include("tests/TransientCheckingpointingTests/TransientLinearBoussinesq.jl") 

    with_mpi() do distribute
        main_transient(distribute,204;restart=false,n_ref_lvls=5,CFL=0.4)
    end     
' > /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
