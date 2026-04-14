using DrWatson
# add Gridap#master
# add GridapDistributed#master

using GridapP4est
using Gridap
using MPI
using PartitionedArrays

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

function perp(vec)
  VectorValue(-vec[2],vec[1])
end

f_scalar(x) = cos(2*π*x[1])*cos(2*π*x[2])
gradperp_f(x) = perp(gradient(f_scalar)(x))


writevtk(Triangulation(models[2]),dir*"/scalar",
cellfields=["f"=>f_scalar, "gradPerp"=>gradperp_f, "div"=>divergence(gradperp_f)],append=false);

ps = [1,2]
ls = LUSolver()
simName ="convergence"

## test convergence
p_convergence_test(ranks,ps,models,H1,dir,f_scalar,gradperp_f,ls,false,true)
plot_convergence_from_saved(dir,simName,["Interp","perp"])

## test div
p_convergence_test(ranks,ps,models,H1,dir,f_scalar,gradperp_f,ls,true,true)
plot_convergence_from_saved(dir,simName,["sum perp","m"])


function H1(
  model::OctreeDistributedDiscreteModel,
  p_fe::Int,dir::String,f_scalar::Function,gradperp_f::Function,ls=LUSolver(),conff=false,return_vtk=false)

  lvl = Int(floor(log2(num_cells(model))))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,4)


  V = FESpace(model, ReferenceFE(raviart_thomas, Float64, p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  P = FESpace(model, ReferenceFE(lagrangian, Float64, p_fe+1); conformity=:H1)
  H = TrialFESpace(P)

  Πperp = interpolate(gradperp_f,U)
  sum(∫( divergence(Πperp) )dΩ)

  fh = interpolate(f_scalar,H)
  gradperp_fh = perp∘(gradient(fh))
  # π0u = interpolate(u_scalar,U)
  # grap_perp_π0u = perp∘(∇(π0u))
  conf = abs(sum.(sum(∫( (perp∘(∇(fh)))  )dΩ )))

  e = f_scalar - fh
  e_f = sqrt(sum(∫( e⋅e )dΩ))

  e_perp = gradperp_fh - Πperp
  e_hdiv = sqrt(sum(∫( e_perp⋅e_perp )dΩ))

  if return_vtk
    writevtk(Ω,dir*"/Hdiv_nref$(lvl)_p$p_fe",
      cellfields=["fh"=>fh,"error_f"=>fh-f_scalar,
              "gradperp"=>Πperp,"gradpherh"=>gradperp_fh,
              "error_gradperp"=>Πperp-gradperp_fh],append=false);
  end


  if conff
    return  conf, false, false
  else
    return e_f, e_hdiv, false
  end
end
