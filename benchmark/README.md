# Benchmarking  

This suite uses `GFlops.jl` and `BenchmarkTools.jl` to compare the performance of computing the grad-grad term for the extrinsic and intrinsic models. 

### Running Locally

To run `test/benchmark.jl` locally, ensure you have installed the packages as follows: 

``` 
 pkg> add GridapGeosciences#master
 pkg> add https://github.com/tamaratambyah/GFlops.jl#remove_remfloat
```
In [my fork of GFlops](https://github.com/tamaratambyah/GFlops.jl/tree/remove_remfloat), I removed the issue with remfloat and output more features from the @gflops macro. 
