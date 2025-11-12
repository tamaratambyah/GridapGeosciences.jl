#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=6
#PBS -l mem=24gb
#PBS -N linear_boussineq
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 6 julia --project=$PBS_O_WORKDIR -e'
    using DrWatson
    include(projectdir("tests/Geophysical/TransientLinearBoussinesq.jl"))

' 
