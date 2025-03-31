using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")


manifold_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
ref_manifold_model = Adaptivity.refine(manifold_model)
ref_ref_manifold_model = Adaptivity.refine(ref_manifold_model)

writevtk(get_ambient_grid(get_grid(ref_manifold_model)),dir*"/ref_grid",append=false)
writevtk(get_ambient_grid(get_grid(ref_ref_manifold_model)),dir*"/ref_ref_grid",append=false)

Ω = Triangulation(manifold_model)
ref_Ω = Triangulation(ref_manifold_model)
ref_ref_Ω = Triangulation(ref_ref_manifold_model)
