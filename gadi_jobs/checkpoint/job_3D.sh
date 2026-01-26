#!/bin/bash
#PBS -P zg98
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l ncpus=208
#PBS -l mem=1024gb
#PBS -N linear_boussineq_amg
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 192 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays 
    include("tests/TransientCheckingpointingTests/TransientLinearBoussinesq.jl") 

    with_mpi() do distribute
        main_transient(distribute,192;restart=false,n_ref_lvls=4,CFL=0.2)
    end     
' > /scratch/$PROJECT/tt4814/LB_3D/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/LB_3D/${PBS_JOBNAME}.err.${PBS_JOBID}
