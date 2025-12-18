using MPI
using PartitionedArrays
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson
using Test

dir = datadir("Nonconfirming")
!isdir(dir) && mkdir(dir)

MPI.Init()
np = MPI.Comm_size(MPI.COMM_WORLD)
ranks = distribute_with_mpi(LinearIndices((np,)))

order = 2
function u_ex(x) ### this exact solution is zero mean + periodic on [0,1]^2
  if x[1] < 0.5
    return x[1]*(0.5-x[1])
  else
    return (x[1]-0.5)*(x[1]-1)
  end
end
f_ex(x) = -1.0*Δ(u_ex)(x)

coarse_model = CartesianDiscreteModel((0,1,0,1),(10,10),isperiodic=(true,true))
dmodel = OctreeDistributedDiscreteModel(ranks,coarse_model)#,2)

ref_coarse_flags=map(partition(get_cell_gids(dmodel.dmodel)), local_views(dmodel) ) do indices,lmodel
    flags=zeros(Cint,length(indices))
    flags.=nothing_flag

    cmap = get_cell_map(get_grid(lmodel))
    ref_points = get_cell_ref_coordinates(lmodel)
    coords = lazy_map(evaluate,cmap,ref_points)

    for (i,xy) in enumerate(coords)
      x = map(x->x[1],xy)
      y = map(x->x[2],xy)
      if any(x .> 0.4 ) && any(x .< 0.6 ) && any(y .> 0.4) && any(y .< 0.6)
          flags[i] = refine_flag
      end
    end
    flags
end
fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

##### Choose a model
model = fmodel
# model = dmodel
# model = coarse_model


#### Triangulate and check compatibility
Ω = Triangulation(model)
degree = 2*(order+1)
dΩ = Measure(Ω,degree)

# check zero mean
sum(∫(u_ex)dΩ  )

# check compatibility
sum(∫( laplacian(u_ex))dΩ  )

################################################################################
####### use zeromean constraint in FE space
V0 = FESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:H1)
# V = Gridap.FESpaces.ZeroMeanFESpace(V0,dΩ)
V = FESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:H1,constraint=:zeromean)
U = TrialFESpace(V)

map(V.spaces) do space
  typeof(space)<: Gridap.FESpaces.FESpaceWithLinearConstraints
end

# Gridap.FESpaces.get_cell_constraints(f::Gridap.FESpaces.FESpaceWithConstantFixed,x::Gridap.FESpaces.Constrained) = Gridap.FESpaces.get_cell_constraints(f.space,x)
# Gridap.FESpaces.get_cell_constraints(f::Gridap.FESpaces.FESpaceWithLinearConstraints,x::Gridap.FESpaces.Constrained) = Gridap.FESpaces.get_cell_constraints(f)
# Gridap.FESpaces.get_cell_isconstrained(f::Gridap.FESpaces.FESpaceWithConstantFixed,x::Gridap.FESpaces.Constrained) = Gridap.FESpaces.get_cell_isconstrained(f.space,x)
# Gridap.FESpaces.get_cell_isconstrained(f::Gridap.FESpaces.FESpaceWithLinearConstraints,x::Gridap.FESpaces.Constrained) = Gridap.FESpaces.get_cell_isconstrained(f)

Gridap.FESpaces.get_cell_constraints(f::Gridap.FESpaces.FESpaceWithConstantFixed) = Gridap.FESpaces.get_cell_constraints(f.space)
Gridap.FESpaces.get_cell_isconstrained(f::Gridap.FESpaces.FESpaceWithConstantFixed) = Gridap.FESpaces.get_cell_isconstrained(f.space)

biform(u,v) = ∫( ∇(u)⋅∇(v)  )dΩ
liform(v) = ∫( v*f_ex )dΩ
op = AffineFEOperator(biform,liform,U,V)
A = get_matrix(op)
b = get_vector(op)

using LinearAlgebra
evals = eigvals(Array(partition(A).item))
# # println(A*(3*ones(size(b))))
sum(b)

uh = solve(LUSolver(),op)

l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
eu = u_ex - uh
eu_l2 = l2(eu)

writevtk(Ω,dir*"/poisson_sol",cellfields=["u"=>u_ex,"uh"=>uh,"eu"=>uh-u_ex],append=false)

################################################################################
####### use lagrange multiplers to enforce zeromean
V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order); conformity=:H1)
U = TrialFESpace(V)

Λ = ConstantFESpace(model) ### broken on trian, must use model
M = TrialFESpace(Λ)

X = MultiFieldFESpace([U,M])
Y = MultiFieldFESpace([V,Λ])

###
poisson_biformX((u,μ),(v,λ)) = ∫( ∇(u)⋅∇(v)  )dΩ  + ∫(v*μ)dΩ + ∫(λ*u)dΩ
poisson_liformY((v,λ)) =  ∫( f_ex*v )dΩ  + ∫(λ*u_ex)dΩ

op = AffineFEOperator(poisson_biformX,poisson_liformY,X,Y)

## hack around the issue https://github.com/tamaratambyah/GridapGeosciences.jl/issues/5
function Gridap.CellData.get_triangulation(a::GridapDistributed.DistributedMultiFieldCellField)
  trians = map(get_triangulation,a.field_fe_fun)
  # @check all(map(t -> t === first(trians), trians))
  return first(trians)
end

uh,μh = solve(LUSolver(),op)

#### Compute errors
e = uh-u_ex
sqrt(sum(∫(e⊙e)dΩ))
writevtk(Ω,dir*"/poisson_sol_lagrange",cellfields=["u"=>u_ex,"uh"=>uh,"eu"=>uh-u_ex],append=false)
