#!/bin/bash
#PBS -P bt62
#PBS -q normalsr
#PBS -l walltime=24:00:00
#PBS -l ncpus=208
#PBS -l mem=1000gb
#PBS -N linear_boussineq_amg
#PBS -l wd

source $HOME/scripts/load-configs-bt62.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 192 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays 
    include("test/TransientCheckingpointingTests/TransientLinearBoussinesq.jl") 

    MPI.Init()
    np = MPI.Comm_size(MPI.COMM_WORLD)
    ranks = distribute_with_mpi(LinearIndices((np,)))

    with_mpi() do distribute
        main_transient(distribute,np;restart=true,n_ref_lvls=4,CFL=0.2)
    end     
' > /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.out.${PBS_JOBID} 2> /scratch/$PROJECT/tt4814/${PBS_JOBNAME}.err.${PBS_JOBID}
