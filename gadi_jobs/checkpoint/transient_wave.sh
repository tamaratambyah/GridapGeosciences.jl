#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=4
#PBS -l mem=16gb
#PBS -N wave_n1
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n $PBS_NCPUS julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays

    MPI.Init()
    nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
    include("tests/TransientCheckingpointingTests/TransientWaveEquation.jl") 

    with_mpi() do distribute
        main_transient(distribute,nprocs;restart=false,n_ref_lvls=1,p_fe=1)
    end

' 
