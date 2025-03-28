using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Helpers
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

using Plots
include("initialise.jl")


for order in [1,2,3,4,5]
  errs = []
  ns = []
  model = ref_ref_ref_ref_ref_model
  while Gridap.Adaptivity.is_child(model,model.parent) && (typeof(model.parent) <: Gridap.Adaptivity.AdaptedDiscreteModel)
    manifold_model_h = ManifoldDiscreteModel(model,cubedsphere)
    nch = sqrt(num_cells(manifold_model_h)/6)
    Ωh = Triangulation(manifold_model_h)
    mh = MetricInfo(metric_func,Ωh)
    dΩgh = Measure(mh,Ωh,order)
    err_h = sum( integrate(1.0,dΩgh)) - 4*π*r^2

    manifold_model_H = ManifoldDiscreteModel(model.parent,cubedsphere)
    ncH = sqrt(num_cells(manifold_model_H)/6)
    ΩH = Triangulation(manifold_model_H)
    mH = MetricInfo(metric_func,ΩH)
    dΩgH = Measure(mH,ΩH,order)
    err_H = sum( integrate(1.0,dΩgH)) - 4*π*r^2



    push!(errs, log2(err_H/err_h) )
    push!(ns,log2(nch))
    model = model.parent
  end

  plot(ns,errs,lw=2,
      xlabel="log2(nc_H)",
      ylabel="log2(err_H/err_h)",
      title = "quad order = $order")
  savefig(plotsdir()*"/order_$order")

end
