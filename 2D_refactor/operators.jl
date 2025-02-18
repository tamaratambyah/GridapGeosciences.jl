using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
using StaticArrays: SVector

include("cube_surface_1_cell_per_panel.jl")
include("structs.jl")
include("refinement.jl")
include("maps.jl")




function metric(θϕ)
  θ,ϕ = θϕ
  TensorValue{2,2}( (cos(ϕ))^2, 0,
                      0,        1)

  # TensorValue{2,2}( 1, 0, 0, 1)
end

"""
surface_gradient -- function input
"""
function surface_gradient(f::Function,g::Function)
  println("function surface grad")
  function _surface_gradient(θϕ)
    surface_gradient(f,g,θϕ)
  end
end

function surface_gradient(f::Function,g::Function,θϕ::Point)
  println("function surface grad point")

  _g=g(θϕ)
  inv(_g)⋅( gradient(f,θϕ) )
end

function surface_gradient(f::Function,g::Function,θϕ::SVector)
  println("function surface grad svector")

  surface_gradient(f,g,Point(θϕ))
end

"""
surface_gradient -- cellfield input
"""
function surface_gradient(a::CellField,g::Function)
  println("cellfield function surface grad")
  _g = CellField(g,get_triangulation(a))
  surface_gradient(a,_g)
end

function surface_gradient(a::CellField,g::CellField)
  println("cellfield cellfield surface grad")

  println(inv(g))
  ( gradient(a) )
  # inv(g)⋅( gradient(a) )
end


function surface_divergence(u::Function,g::Function)
  function _surface_divergence(θϕ)
     _g=g(θϕ)
     function f(θϕ)
        sqrt(det(_g))*(u(θϕ))
     end
     1.0/sqrt(det(_g))*(∇⋅(f))(θϕ)
  end
end


function surface_laplacian(f::Function,g::Function)
  function _surface_laplacian(θϕ)
    surface_divergence(surface_gradient(f,g),g)(θϕ)
  end
end

# metric(θϕ)
u(θϕ) = 3*θϕ[1] + 2*θϕ[2]^2
# w(θϕ) = VectorValue( 3*θϕ[1],2*θϕ[2])
# θϕ = VectorValue(π,π/4)
# evaluate(surface_gradient(u,metric),θϕ)
# evaluate(surface_divergence(w),θϕ)
# evaluate(surface_laplacian(u),θϕ)

f(θϕ) = surface_laplacian(u,metric)(θϕ)

# a(R) = R/sqrt(3)
# R(a) = a*sqrt(3)

# function A_other(θϕ,a)
#   θ,ϕ = θϕ
#   (R(a)*cos(θ)*cos(ϕ)/a)*( TensorValue{2,2}( cos(θ), 0,
#                   -sin(θ)*sin(ϕ), cos(ϕ)) )
# end

# function A_top(θϕ,a)
#   θ,ϕ = θϕ
#   (R(a)*sin(ϕ)/a)*( TensorValue{2,2}( cos(θ), sin(θ),
#                   -sin(θ)*sin(ϕ), sin(ϕ)*cos(θ) ) )

# end

# function A_bottom(θϕ,a)
#   θ,ϕ = θϕ
#   (R(a)*sin(ϕ)/a)*( TensorValue{2,2}( -cos(θ), sin(θ),
#                   sin(θ)*sin(ϕ), sin(ϕ)*cos(θ) ) )

# end

# function metric(θϕ)
#   X,Y,Z =  θϕ2xyz(θϕ)
#   a=1
#   println(Z)
#   if Z == 1
#     A = A_top(θϕ,a)
#   elseif Z == -1
#     A = A_bottom(θϕ,a)
#   else
#     A = A_other(θϕ,a)
#   end
#   g = transpose(A)⋅A
#   println(A)
#   println(g)
#   println(inv(g))
#   g
# end


t = evaluate(metric,get_cell_points(Ω))






dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
CSgrid = CubedSphereGrid(cube_grid,map_cube_to_latlon)
model = ManifoldDiscreteModel(CSgrid,topo,face_labels)
CSmodelh = Gridap.Adaptivity.refine(model)
CSmodel = Gridap.Adaptivity.refine(CSmodelh)
# model = CartesianDiscreteModel((0,2π,-π/2,π/2), (1,1))

writevtk(model,dir*"/model",append=false)
writevtk(CSmodelh,dir*"/CSmodel",append=false)


p = 1
degree = 2*p+1
V = FESpace(model,ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)



a(u,v) = ∫( surface_gradient(u,metric)⋅surface_gradient(v,metric) )dΩ
_a(u,v) = ∫( gradient(u)⋅gradient(v) )dΩ
b(v)   = ∫( f*v )dΩ

# res(u,v) = a(u,v) - b(v)
# jac(u,du,v) = a(du,v)

# op = AffineFEOperator(a,b,U,V)
# uh = solve(op)

uh = zero(U)
du = get_trial_fe_basis(U)
v = get_fe_basis(V)
assem = SparseMatrixAssembler(U,V)

matdata = collect_cell_matrix(U,V,a(du,v))
A = assemble_matrix(assem, matdata)

_matdata = collect_cell_matrix(U,V,_a(du,v))
_A = assemble_matrix(assem, _matdata)

vecdata = collect_cell_vector(V,b(v))
rhs = assemble_vector(assem, vecdata)

uh = A\rhs

dc = a(du,v)

keys = get_domains(dc)
w = []

for strian in get_domains(dc)
  scell_mat = get_contribution(dc,strian)
  cell_mat, trian = move_contributions(scell_mat,strian)
  push!(w,cell_mat)
end

W = collect(w[1])



# cell_reffe =  Fill(ReferenceFE(lagrangian,Float64,p),num_cells(model))
# conf = H1Conformity()


# cell_fe = Gridap.FESpaces.CellFE(model,cell_reffe,conf)
# trian = Triangulation(model)
# labels = get_face_labeling(model)
# dirichlet_tags=Int[]
# dirichlet_masks=nothing
# constraint=nothing
# vector_type=nothing

# @assert num_cells(cell_fe) == num_cells(model) """\n
# The number of cells provided in the `cell_fe` argument ($(cell_fe.num_cells) cells)
# does not match the number of cells ($(num_cells(model)) cells) in the provided DiscreteModel.
# """
# _vector_type = _get_vector_type(vector_type,cell_fe,trian)
# F = _ConformingFESpace(
#     _vector_type,
#     model,
#     labels,
#     cell_fe,
#     dirichlet_tags,
#     dirichlet_masks,
#     trian)
# V = _add_constraint(F,cell_fe.max_order,constraint)
# V









# solver = NLSolver(LUSolver())
# uh = solve(solver,op)

# e=u-uh
# sqrt(sum(∫( e*e  )dΩ))
