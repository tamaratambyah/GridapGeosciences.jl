#!/bin/bash

for file in gadi_jobs/convergence*.sh; do
  # qsub $file
  echo $file
done