using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.Helpers
using Gridap.TensorValues
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

include("initialise.jl")
include("maps/metric_maps.jl")
include("metric/surface_metric.jl")
include("surface_quadrature.jl")

# _model = cube_model_3D
_model = ref_ref_model
manifold_grid = CubedSphereGrid(_model)
face_labels = get_face_labeling(_model)
order = 4



struct ManifoldDiscreteModel{Dc,Dp,Dp_topo,A<:Grid{Dc,Dp},B<:GridTopology{Dc,Dp_topo}} <: DiscreteModel{Dc,Dp}
  manifold_grid:: A
  grid_topology:: B
  face_labeling::FaceLabeling
end

function ManifoldDiscreteModel(manifold_grid::Grid{Dc,Dp},topo::GridTopology{Dc,Dp_topo},labels) where {Dc,Dp,Dp_topo}
  A = typeof(manifold_grid)
  B = typeof(topo)
  ManifoldDiscreteModel{Dc,Dp,Dp_topo,A,B}(manifold_grid,topo,labels)
end

Gridap.Geometry.get_grid(model::ManifoldDiscreteModel) = model.manifold_grid

Gridap.Geometry.get_grid_topology(model::ManifoldDiscreteModel) = model.grid_topology

Gridap.Geometry.get_face_labeling(model::ManifoldDiscreteModel) = model.face_labeling


manifold_model = ManifoldDiscreteModel(manifold_grid,get_grid_topology(_model),get_face_labeling(_model))
Ω = Triangulation(manifold_model)
m = MetricInfo(metric_func,Ω)
dΩg = Measure(m,Ω,order)
sum( integrate(1.0,dΩg))
4*π*r^2
