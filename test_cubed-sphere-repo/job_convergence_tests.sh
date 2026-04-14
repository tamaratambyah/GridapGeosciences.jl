#!/bin/bash
#PBS -P zg98
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l ncpus=4
#PBS -l mem=16gb
#PBS -N transient_advection_supg
#PBS -l wd

source $HOME/scripts/load-configs-zg98.sh
source $HOME/scripts/load-intel.sh

mpiexec -n 4 julia --project=$PBS_O_WORKDIR -e'
  using MPI
  using PartitionedArrays
  # include("test/Laplace/LaplaceBeltrami.jl")
  # include("test/Laplace/Helmholtz.jl")
  # include("test/Laplace/MixedHelmholtz.jl")
  # include("test/Geophysical/WaveEquation.jl")
  # include("test/Geophysical/ShallowWater.jl")
  # include("test/Geophysical/TransientShallowWater.jl")

  # include("test/Advection/AdvectionSUPG.jl")
  # include("test/Advection/AdvectionDGUpwinding.jl")
  # include("test/Advection/TransientAdvectionDGUpwinding.jl")
  include("test/Advection/TransientAdvectionSUPG.jl")

    with_mpi() do distribute
        # LaplaceBeltrami.main(distribute,4;octree=true)
        # Helmholtz.main(distribute,4;octree=true)
        # MixedHelmholtz.main(distribute,4,true)
        # WaveEquation.main(distribute,4;octree=true)
        # ShallowWater.main(distribute,4;octree=true)
        # TransientShallowWater.main(distribute,4;octree=true)
        # AdvectionSUPG.main(distribute,4;octree=true)
        # AdvectionDGUpwinding.main(distribute,4;octree=true)
        TransientAdvectionSUPG.main(distribute,4;octree=true)
        # TransientAdvectionDGUpwinding.main(distribute,4;octree=true)
    end

' 
