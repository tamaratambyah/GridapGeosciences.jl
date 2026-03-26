#!/bin/bash
#PBS -P bt62
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=192
#PBS -l mem=768gb
#PBS -N sw_3D
#PBS -l wd

source $HOME/scripts/load-configs.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 192 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/TransientCheckingpointingTests/TransientShallowWater_3D.jl") 

    with_mpi() do distribute
        main_transient(distribute,192;restart=true,n_ref_lvls=6,p_fe=0,CFL=0.1)
        # main_visualise(distribute,192;n_ref_lvls=6,p_fe=0)
    end

' > /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
