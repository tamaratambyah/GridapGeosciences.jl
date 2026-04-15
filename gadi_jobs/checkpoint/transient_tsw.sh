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
    include("test_cubed-sphere-repo/TransientCheckingpointingTests/TSWGalewsky.jl") 

    MPI.Init()
    np = MPI.Comm_size(MPI.COMM_WORLD)
    ranks = distribute_with_mpi(LinearIndices((np,)))

    with_mpi() do distribute
        main_transient(distribute,np;restart=false,n_ref_lvls=6,p_fe=1,CFL=0.1)
        # main_visualise(distribute,np;n_ref_lvls=6,p_fe=1)
    end

' > /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
