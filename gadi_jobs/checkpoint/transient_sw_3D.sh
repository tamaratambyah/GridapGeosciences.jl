#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=384
#PBS -l mem=1536gb
#PBS -N sw_3D
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 384 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/TransientCheckingpointingTests/TransientShallowWater_3D.jl") 

    with_mpi() do distribute
        main_transient(distribute,384;restart=true,n_ref_lvls=6,p_fe=0,CFL=0.1)
        # main_visualise(distribute,384;n_ref_lvls=6,p_fe=0)
    end

' > /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
