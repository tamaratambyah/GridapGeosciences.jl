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


function adapt_model(ranks,model; refine=true)
    cell_partition=get_cell_gids(model.octree_model.dmodel)
    ref_coarse_flags=map(ranks,partition(cell_partition)) do rank,indices
        flags=zeros(Cint,length(indices))
        if refine
          flags.=refine_flag
        else
          flags.=nothing_flag
        end
    end
    # Gridap.Adaptivity.adapt(model,ref_coarse_flags)
    GridapP4est.adapt(model,ref_coarse_flags)
end

n_refs = 2 # number of refinements

# coarsest model
model0 = CubedSphereDiscreteModel(ranks,n_refs;adaptive=true)

model1,glue1 = adapt_model(ranks, model0; refine=false)

writevtk(Triangulation(model1),joinpath(datadir("models"),"cubed_sphere_amr_step_0"),append=false)

model2,glue2 = adapt_model(ranks, model1)
writevtk(Triangulation(model2),joinpath(datadir("models"),"cubed_sphere_amr_step_1"),append=false)

models = [model2,model1]
glues = [glue2,glue1]

######

T = Int64
coarsest_model = CubedSphereDiscreteModel(ranks,n_refs;adaptive=true)
num_procs_x_level = [(1,1),(1,1)]
_num_levels  = length(num_procs_x_level)
level_parts = Vector{Union{typeof(ranks),GridapDistributed.MPIVoidVector{T}}}(undef,_num_levels)
meshes      = Vector{ModelHierarchyLevel}(undef,_num_levels)

level_parts[_num_levels] = get_parts(coarsest_model)
meshes[_num_levels] = ModelHierarchyLevel(_num_levels,coarsest_model,nothing,nothing,nothing)

i = 1

level_parts[i]     = level_parts[i+1]
model_ref,ref_glue = adapt_model(ranks,coarsest_model) #Gridap.Adaptivity.refine(modelH)
model_red,red_glue = nothing,nothing
meshes[i] = ModelHierarchyLevel(i,model_ref,ref_glue,model_red,red_glue)

mh = HierarchicalArray(meshes,level_parts)

model     = get_model_before_redist(model,lev)


_mh1 = map(linear_indices(mh),mh) do lev, mhl
  if lev == num_levels(mh)
    return mhl
  end

  if GridapDistributed.i_am_in(get_level_parts(mh,lev+1))
    model     = get_model_before_redist(mh,lev)
    println(model)
    parent    = get_model(mh,lev+1)
    ref_glue  = mhl.ref_glue
    model_ref = GridapDistributed.DistributedAdaptedDiscreteModel(model,parent,ref_glue)
    # model_ref,ref_glue = adapt_model(ranks,parent)
  else
    model_ref = nothing
  end
  return ModelHierarchyLevel(lev,model_ref,mhl.ref_glue,mhl.model_red,mhl.red_glue)
end



# function CubedSphereModelHierarchy(models,glues)
#   nlevs = length(models)
#   @check length(models) == length(glues) "Incorrect dimensions"
#   # for lev in 1:nlevs-1
#   #   @check num_cells(models[lev]) > num_cells(models[lev+1]) "Incorrect hierarchy of models."
#   # end
#   ranks = get_parts(models[1])
#   @check all(m -> length(get_parts(m)) === length(ranks), models) "Models have different communicators."

#   level_parts = fill(ranks,nlevs)
#   meshes = Vector{ModelHierarchyLevel}(undef,nlevs)
#   for lev in 1:nlevs-1
#     model = models[lev]
#     glue  = glues[lev] #Gridap.Adaptivity.get_adaptivity_glue(models[lev])

#     meshes[lev] = ModelHierarchyLevel(lev,model,glue,nothing,nothing)
#   end
#   meshes[nlevs] = ModelHierarchyLevel(nlevs,models[nlevs],nothing,nothing,nothing)
#   return HierarchicalArray(meshes,level_parts)
# end


# mh = CubedSphereModelHierarchy(models,glues)

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



model = get_model(_mh1,1)
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

reffe_u = ReferenceFE(raviart_thomas,Float64,p)
Vs = TestFESpace(_mh1,reffe_u,conformity=:Hdiv)
Us = TrialFESpace(Vs)

reffe_p = ReferenceFE(lagrangian,Float64,p)
Qs = TestFESpace(_mh1,reffe_p;conformity=:L2)
Ps = TrialFESpace(Qs)

tests = MultiFieldFESpace([Us,Ps])
trials = MultiFieldFESpace([Vs,Qs])


biform((u,p),(v,q),dΩ) = ( ∫( u⋅v - (∇⋅v)*p )dΩ
                         + ∫( p*q + (∇⋅u)*q )dΩ )
liform((v,q),dΩ) = ∫( f_u⋅v  + f_p⋅q )dΩ

a(u,v) = biform(u,v,dΩ)
l(v) = liform(v,dΩ)
op = AffineFEOperator(a,l,get_fe_space(trials,1),get_fe_space(tests,1))
A = get_matrix(op)
b = get_vector(op)

biforms = map(mhl -> get_bilinear_form(mhl,biform,degree),_mh1)

# smoothers = get_patch_smoothers(mh,tests,biform,degree)
smoothers = get_jacobi_smoothers(_mh1)

restrictions, prolongations = setup_transfer_operators(
  trials, degree; mode=:residual, solver=GridapSolvers.LinearSolvers.IS_ConjugateGradientSolver(;reltol=1.e-6)
)


gmg = GMGLinearSolver(
  _mh1,trials,tests,biforms,
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
