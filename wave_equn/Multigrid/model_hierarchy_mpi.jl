using GridapGeosciences
using GridapP4est
using PartitionedArrays
using MPI
using Gridap
using GridapDistributed
using DrWatson
using FillArrays
using GridapSolvers
using Gridap.Helpers: @check
import GridapSolvers.MultilevelTools: ModelHierarchyLevel, HierarchicalArray


nprocs = (1,1)
ranks = with_mpi() do distribute
  distribute(LinearIndices((prod(nprocs),)))
end


function adapt_model(ranks,model)
    cell_partition=get_cell_gids(model.octree_model.dmodel)
    ref_coarse_flags=map(ranks,partition(cell_partition)) do rank,indices
        flags=zeros(Cint,length(indices))
          flags.=refine_flag
    end
    Gridap.Adaptivity.adapt(model,ref_coarse_flags)
    # GridapP4est.adapt(model,ref_coarse_flags)
end

n_refs = 1 # number of refinements

T = Int64
coarsest_model = CubedSphereDiscreteModel(ranks,n_refs;adaptive=true)
num_procs_x_level = [(1,1),(1,1)]

_num_levels  = length(num_procs_x_level)
level_parts = Vector{Union{typeof(ranks),GridapDistributed.MPIVoidVector{T}}}(undef,_num_levels)
meshes      = Vector{ModelHierarchyLevel}(undef,_num_levels)

level_parts[_num_levels] = get_parts(coarsest_model)
meshes[_num_levels] = ModelHierarchyLevel(_num_levels,coarsest_model,nothing,nothing,nothing)

i = _num_levels-1

level_parts[i]     = level_parts[i+1]
model_ref,ref_glue = adapt_model(ranks,coarsest_model)   #Gridap.Adaptivity.refine(modelH)
model_red,red_glue = nothing,nothing
meshes[i] = ModelHierarchyLevel(i,model_ref,ref_glue,model_red,red_glue)

mh = HierarchicalArray(meshes,level_parts)


import GridapGeosciences: ForestOfOctreesCubedSphereDiscreteModel
import Gridap.Adaptivity: AdaptivityGlue, AdaptedDiscreteModel
import GridapDistributed: DistributedDiscreteModel
struct ForestAdaptedDiscreteModel{Dc,Dp,A, B,C<:AdaptivityGlue} <: DiscreteModel{Dc,Dp}
  model::A
  parent::B
  glue::C
  function ForestAdaptedDiscreteModel(model,parent,glue)
    # @check !isa(model,AdaptedDiscreteModel)
    A = typeof(model)
    B = typeof(parent)
    C = typeof(glue)
    Dc = 2
    Dp = 3
    return new{Dc,Dp,A,B,C}(model,parent,glue)
  end
end


# DiscreteModel API
Gridap.Geometry.get_grid(model::ForestAdaptedDiscreteModel)          = get_grid(model.model)
Gridap.Geometry.get_grid_topology(model::ForestAdaptedDiscreteModel) = Gridap.Geometry.get_grid_topology(model.model)
Gridap.Geometry.get_face_labeling(model::ForestAdaptedDiscreteModel) = get_face_labeling(model.model)

# Other getters
Gridap.Adaptivity.get_model(model::ForestAdaptedDiscreteModel)  = model.model
Gridap.Adaptivity.get_parent(model::ForestAdaptedDiscreteModel{Dc,Dp,A,<:ForestAdaptedDiscreteModel}) where {Dc,Dp,A} = get_model(model.parent)
Gridap.Adaptivity.get_parent(model::ForestAdaptedDiscreteModel{Dc,Dp,A,B}) where {Dc,Dp,A,B} = model.parent
Gridap.Adaptivity.get_adaptivity_glue(model::ForestAdaptedDiscreteModel) = model.glue

function is_child(m1::ForestAdaptedDiscreteModel,m2)
  return get_parent(m1) === m2 # m1 = refine(m2)
end

function is_child(m1::ForestAdaptedDiscreteModel,m2::ForestAdaptedDiscreteModel)
  return get_parent(m1) === get_model(m2) # m1 = refine(m2)
end

is_child(m1,m2::ForestAdaptedDiscreteModel) = false

is_related(m1,m2) = is_child(m1,m2) || is_child(m2,m1)


function GridapDistributed.DistributedAdaptedDiscreteModel(
  model  :: DistributedDiscreteModel,
  parent :: DistributedDiscreteModel,
  glue   :: AbstractArray{<:AdaptivityGlue};
)
  models = map(local_views(model),local_views(parent),glue) do model, parent, glue
    ForestAdaptedDiscreteModel(model,parent,glue)
  end
  gids = get_cell_gids(model)
  metadata = hasproperty(model,:metadata) ? model.metadata : nothing
  return GridapDistributed.GenericDistributedDiscreteModel(models,gids;metadata)

end


mh = GridapSolvers.MultilevelTools.convert_to_adapted_models(mh)
lev = 1
mhl = ModelHierarchyLevel(i,model_ref,ref_glue,model_red,red_glue)
model     = get_model_before_redist(mh,lev)
parent    = get_model(mh,lev+1)
ref_glue  = mhl.ref_glue




function get_jacobi_smoothers(mh)
  nlevs = num_levels(mh)
  smoothers = Fill(RichardsonSmoother(JacobiLinearSolver(),10,2.0/3.0),nlevs-1)
  level_parts = view(get_level_parts(mh),1:nlevs-1)
  return HierarchicalArray(smoothers,level_parts)
end

function get_patch_smoothers(mh,tests,biform,qdegree,ω=0.2,niter=10,local_solver=LUSolver(),is_nonlinear=false,weighted=false)
  patch_decompositions = PatchDecomposition(mh)
  patch_spaces = PatchFESpace(tests,patch_decompositions)
  nlevs = num_levels(mh)
  smoothers = map(view(tests,1:nlevs-1),patch_decompositions,patch_spaces) do tests, PD, Ph
    Vh = get_fe_space(tests)
    Ω  = Triangulation(PD)
    dΩ = Measure(Ω,qdegree)
    ap = (u,v) -> biform(u,v,dΩ)
    patch_smoother = GridapSolvers.PatchBasedSmoothers.PatchBasedLinearSolver(ap,Ph,Vh;local_solver,is_nonlinear,weighted)
    return RichardsonSmoother(patch_smoother,niter,ω) # 10 iterations, ω = 0.2
  end
  return smoothers
end

function get_bilinear_form(mh_lev,biform,qdegree)
  model = get_model(mh_lev)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,qdegree)
  return (u,v) -> biform(u,v,dΩ)
end

u_exact(x) = VectorValue(x[1]*(1-x[1]), x[2]*(1-x[2]),0.0)
p_exact(x) = 1.0 + 0.001*x[1]*(1-x[1]) + 0.001*x[2]*(1-x[2])
f_u(x) =  u_exact(x) + ∇(p_exact)(x)
f_p(x) = p_exact(x) + (∇⋅u_exact)(x)


p = 1
degree = 6

import Gridap.Geometry: Grid, GridTopology
Gridap.Geometry.num_point_dims(stuff::Grid) = 3
Gridap.Geometry.num_point_dims(model::DiscreteModel) = 3
Gridap.Geometry.num_point_dims(::GridTopology) =3

model = get_model(mh,1)
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

reffe_u = ReferenceFE(raviart_thomas,Float64,p)
Vs = TestFESpace(mh,reffe_u,conformity=:Hdiv)
Us = TrialFESpace(Vs)

reffe_p = ReferenceFE(lagrangian,Float64,p)
Qs = TestFESpace(mh,reffe_p;conformity=:L2)
Ps = TrialFESpace(Qs)

# tests = MultiFieldFESpace([Us,Ps])
# trials = MultiFieldFESpace([Vs,Qs])
_biform(p,q,dΩ) =  ∫( p*q  )dΩ
_liform(q,dΩ) = ∫( 0.0*q )dΩ


# biform((u,p),(v,q),dΩ) = ( ∫( u⋅v - (∇⋅v)*p )dΩ
#                          + ∫( p*q + (∇⋅u)*q )dΩ )
# liform((v,q),dΩ) = ∫( f_u⋅v  + f_p⋅q )dΩ

a(u,v) = _biform(u,v,dΩ)
l(v) = _liform(v,dΩ)
op = AffineFEOperator(a,l,get_fe_space(Ps,1),get_fe_space(Qs,1))
A = get_matrix(op)
b = get_vector(op)

biforms = map(mhl -> get_bilinear_form(mhl,biform,degree),mh)

# smoothers = get_patch_smoothers(mh,tests,biform,degree)
smoothers = get_jacobi_smoothers(mh)

restrictions, prolongations = setup_transfer_operators(
  trials, degree; mode=:residual, solver=GridapSolvers.LinearSolvers.IS_ConjugateGradientSolver(;reltol=1.e-6)
)


gmg = GMGLinearSolver(
  mh,trials,tests,biforms,
  prolongations,restrictions,
  pre_smoothers=smoothers,
  post_smoothers=smoothers,
  coarsest_solver=LUSolver(),
  maxiter=3,mode=:preconditioner,
  verbose=true
)
gmg.log.depth = 2


# smatrices, A, b = compute_hierarchy_matrices(trials,tests,biform,liform,degree)

solver = GMRESSolver(20;Pl=gmg,maxiter=5,atol=1e-14,rtol=1.e-14,verbose=true)
ns = numerical_setup(symbolic_setup(solver,A),A)

# Solve
x = pfill(0.0,partition(axes(A,2)))
solve!(x,ns,b)


# Error
Uh = get_fe_space(trials,1)
uh, ph = FEFunction(Uh,x)

l2(u_exact-uh,dΩ)
l2(p_exact-ph,dΩ)

writevtk(Ω,datadir("wave"),cellfields = ["uh"=>uh, "ph"=>ph,
                "u"=>CellField(u_exact,Ω), "p"=>CellField(p_exact,Ω)],append=false)
