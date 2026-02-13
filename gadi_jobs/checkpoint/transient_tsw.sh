#!/bin/bash
#PBS -P bt62
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=192
#PBS -l mem=768gb
#PBS -N tsw_galewsky
#PBS -l wd

source $HOME/scripts/load-configs.sh
source $HOME/scripts/load-intel.sh

mpiexec -n $PBS_NCPUS julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/TransientCheckingpointingTests/TSWGalewsky.jl") 

    with_mpi() do distribute
        main_transient(distribute,192;restart=false,n_ref_lvls=5,p_fe=1)
    end

' > /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
