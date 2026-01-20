#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=48
#PBS -l mem=192gb
#PBS -N sw_n3
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n $PBS_NCPUS julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/TransientCheckingpointingTests/TransientShallowWater.jl") 

    with_mpi() do distribute
        main_transient(distribute,48;restart=false,n_ref_lvls=3,p_fe=2)
    end

' > /scratch/$PROJECT/tt4814/SW_W2_convergence_p2/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/SW_W2_convergence_p2/${PBS_JOBNAME}.err.${PBS_JOBID}
