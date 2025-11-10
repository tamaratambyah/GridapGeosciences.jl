#!/bin/bash
files=(
    gadi_jobs/convergence_laplace.sh
    gadi_jobs/convergence_wave.sh
)

# for file in gadi_jobs/convergence*.sh; do
for file in "${files[@]}"; do
  qsub $file
  # echo $file
done