using GridapP4est
# add Gridap#master
using Gridap
using MPI
using PartitionedArrays
using DrWatson
using GridapDistributed

include("refinement_helpers.jl")
include("../convergence_tools.jl")


dir = datadir("derham_trig")
!isdir(dir) && mkdir(dir)

n_ref_lvls = 3

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
coarse_model = CartesianDiscreteModel((0,1,0,1),(4,4),isperiodic=(true,true))
dmodel = OctreeDistributedDiscreteModel(ranks, coarse_model)
ref_coarse_flags = initial_unbalance(dmodel)
fmodel, =Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

models = Vector{OctreeDistributedDiscreteModel}(undef,n_ref_lvls+1)
models[end] = fmodel
for (i,n) in enumerate(n_ref_lvls:-1:1)
  fmodel = refine_all(fmodel)
  models[n] = fmodel
end

# # plot all models
# for i in 1:length(models)
#   ref = length(models)-i
#   writevtk(Triangulation(models[i]),dir*"/model_$ref",append=false);
# end

u_vector(x) = VectorValue(sin(2*π*x[1])*cos(2*π*x[2]), -sin(2*π*x[2])*cos(2*π*x[1]))
div_u(x) = (∇⋅u_vector)(x)

writevtk(Triangulation(models[2]),dir*"/velocity",
cellfields=["u"=>u_vector, "div_u"=>div_u],append=false);

ps = [1]
ls = LUSolver()
simName ="convergence"

## test convergence
p_convergence_test(ranks,ps,models,Hdiv,dir,u_vector,div_u,ls,false,true)
plot_convergence_from_saved(dir,simName,["Interp","div"])

## test div
p_convergence_test(ranks,ps,models,Hdiv,dir,u_vector,div_u,ls,true,true)
plot_convergence_from_saved(dir,simName,["sum div","m"])


function Hdiv(
  model::OctreeDistributedDiscreteModel,
  p_fe::Int,dir::String,u_vector::Function,div_u::Function,ls=LUSolver(),conff=false,return_vtk=false)

  lvl = Int(floor(log2(num_cells(model))))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,4)


  V = FESpace(model, ReferenceFE(raviart_thomas, Float64, p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  W = FESpace(model, ReferenceFE(lagrangian, Float64, p_fe); conformity=:L2)
  R = TrialFESpace(W)

  Πdiv_u = interpolate(div_u,R)
  sum(∫( div_u )dΩ)
  sum(∫( Πdiv_u  )dΩ )

  # u_h = interpolate(u_vector,U)
  a(u,v) = ∫(u⋅v)dΩ
  l(v) = ∫( u_vector⋅v )dΩ
  op = AffineFEOperator(a,l,U,V)
  u_h = solve(ls,op)

  div_uh = divergence(u_h)
  # _div_uh = interpolate(div_uh,R)
  conf = abs(sum(∫( div_uh )dΩ ))
  # sum(∫( _div_uh )dΩ )

  e = u_vector -u_h
  e_u = sqrt(sum(∫( e⋅e )dΩ))

  e_div = div_uh - Πdiv_u
  e_hdiv = sqrt(sum(∫( e_div⋅e_div )dΩ))

  if return_vtk
    writevtk(Ω,dir*"/Hdiv_nref$(lvl)_p$p_fe",
      cellfields=["uh"=>u_h,"error_u"=>u_h-u_vector,
              "div_u"=>Πdiv_u,"div_uh"=>div_uh,
              "error_div"=>Πdiv_u-div_uh],append=false);
  end


  if conff
    return  conf, false, false
  else
    return e_u, e_hdiv, false
  end
end
