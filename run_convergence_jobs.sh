#!/bin/bash
files=(
    gadi_jobs/convergence/convergence_laplace.sh
    gadi_jobs/convergence/convergence_wave.sh
    gadi_jobs/convergence/convergence_linear_shallow_water.sh
    gadi_jobs/convergence/convergence_boussinesq.sh
)

# for file in gadi_jobs/convergence*.sh; do
for file in "${files[@]}"; do
  qsub $file
  # echo $file
done