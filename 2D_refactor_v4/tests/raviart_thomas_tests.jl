using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays


include("../src/initialise.jl")



######### cubed sphere model
coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
model = Adaptivity.refine(coarse_model)
# ref_model = Adaptivity.refine(Adaptivity.refine(model))


V = TestFESpace(coarse_model,ReferenceFE(raviart_thomas,Float64,order); conformity=:Hdiv)

function Gridap.FESpaces.get_cell_dof_basis(
  model::ManifoldDiscreteModel,
  cell_reffe::AbstractArray{<:GenericRefFE{<:DivConforming}},
  ::DivConformity,
  sign_flip = FESpaces.get_sign_flip(model, cell_reffe)
)

  println("my func")

  Dc = num_dims(model)
  Dp_topo = num_point_dims(get_parametric_grid(get_grid(model)))

  cell_map  = get_cell_map(get_grid(model))
  Jt  = lazy_map(Broadcasting(∇),cell_map)
  x   = lazy_map(get_nodes,lazy_map(get_dof_basis,cell_reffe))
  Jtx = lazy_map(evaluate,Jt,x)
  k = FESpaces.TransformRTDofBasis{Dc,Dp_topo}()
  lazy_map(k,cell_reffe,Jtx,sign_flip)
end


# model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,2)))

reffes = get_cell_reffe(model)
topo = get_grid_topology(model)
cell_reffe = expand_cell_data(reffes,get_cell_type(topo))

cell_fe = FESpaces.CellFE(model,cell_reffe,Conformity(testitem(cell_reffe)))
FESpaces._get_vector_type(nothing,cell_fe,Triangulation(model))
