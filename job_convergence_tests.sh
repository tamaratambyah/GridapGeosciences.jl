#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l ncpus=4
#PBS -l mem=16gb
#PBS -N convergence
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 4 julia --project=$PBS_O_WORKDIR -e'
  using MPI
  using PartitionedArrays
  include("tests/Laplace/LaplaceBeltrami.jl")
  include("tests/Laplace/Helmholtz.jl")
  include("tests/Laplace/MixedHelmholtz.jl")
  include("tests/Geophysical/WaveEquation.jl")
  include("tests/Geophysical/ShallowWater.jl")
  include("tests/Advection/AdvectionSUPG.jl")
  include("tests/Advection/AdvectionDGUpwinding.jl")
  include("tests/Advection/TransientAdvectionSUPG.jl")

    with_mpi() do distribute
        # LaplaceBeltrami.main(distribute,4;octree=true)
        # Helmholtz.main(distribute,4;octree=true)
        # MixedHelmholtz.main(distribute,4,true)
        # WaveEquation.main(distribute,4;octree=true)
        # ShallowWater.main(distribute,4;octree=true)
        # AdvectionSUPG.main(distribute,4;octree=true)
        # AdvectionDGUpwinding.main(distribute,1)
        TransientAdvectionSUPG.main(distribute,4)
    end

' 
