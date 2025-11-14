#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=6
#PBS -l mem=64gb
#PBS -N linear_boussineq
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 6 julia --project=$PBS_O_WORKDIR -e'
    using MPI
    using PartitionedArrays
    include("tests/Geophysical/TransientLinearBoussinesq.jl")
    
    with_mpi() do distribute
       main_transient(distribute,6;n_ref_lvls=3,p_fe=1,CFL=0.2,return_vtk=true)
    end      
' 
