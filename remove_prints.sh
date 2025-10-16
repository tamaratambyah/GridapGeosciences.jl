#!/bin/bash
#file=src/ODEs/DAE.jl
#sed -i /println/s/^/#/ $file

for file in src/*/*.jl; do
  sed -i /println/s/^/#/ $file
done