using GridapP4est
using Gridap
using MPI
using PartitionedArrays
using DrWatson
using GridapDistributed

include("refinement_helpers.jl")
include("../convergence_tools.jl")

# perp in 2D cartesian space
function perp(vec)
  VectorValue(-vec[2],vec[1])
end

dir = datadir("derham_trig")
!isdir(dir) && mkdir(dir)

# u_scalar(x) = x[1]*(1-x[1])
# u_vector(x) = VectorValue(x[1]*(1-x[1]),x[2]*(1-x[2]))

u_scalar(x) = sin(2*π*x[1])
u_vector(x) = VectorValue(cos(2*π*x[1]), 0.0)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

model0 = CartesianDiscreteModel((0,1,0,1),(4,4),isperiodic=(true,true))
model1 = CartesianDiscreteModel((0,1,0,1),(8,8),isperiodic=(true,true))
model2 = CartesianDiscreteModel((0,1,0,1),(16,16),isperiodic=(true,true))
model3 = CartesianDiscreteModel((0,1,0,1),(32,32),isperiodic=(true,true))
models = [model3, model2, model1, model0]

ps = [1,2]
ls = LUSolver()
p_convergence_test(ranks,ps,models,derham_complex,dir,u_scalar,u_vector,ls,true)

simName ="convergence"
plot_convergence_from_saved(dir,simName,["H1","Hdiv"])

function derham_complex(
  coarse_model::CartesianDiscreteModel,
  p_fe::Int,dir::String,u_scalar::Function,u_vector::Function,ls=LUSolver(),return_vtk=false)

  lvl = nref(num_cells(coarse_model))
  i_am_main(ranks) && println("nref = $lvl; p_fe = $p_fe")


  dmodel = OctreeDistributedDiscreteModel(ranks, coarse_model)
  ref_coarse_flags = middle_refinement(dmodel)
  fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

  model = fmodel

  V = FESpace(model,ReferenceFE(lagrangian, Float64, p_fe+1); conformity=:H1)
  U = TrialFESpace(V)

  R = FESpace(model, ReferenceFE(raviart_thomas, Float64, p_fe); conformity=:HDiv)
  T = TrialFESpace(R)

  D = FESpace(model, ReferenceFE(lagrangian, Float64, p_fe); conformity=:L2)
  E = TrialFESpace(D)

  degree = 6*(p_fe+1)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  # interpolate scalar into H1, take perp, interpolate perp into RT
  π0u = interpolate(u_scalar,U)
  grap_perp_π0u = perp∘(∇(π0u))
  π1u_grad_perp_u = interpolate(perp∘∇(u_scalar),R)

  e = π1u_grad_perp_u-grap_perp_π0u
  e_h1 = sqrt(sum(∫(e⋅e)dΩ))


  a(u,v) = ∫( u⋅v )dΩ
  l(v) = ∫( u_vector⋅v )dΩ
  op = AffineFEOperator(a,l,T,R)
  π1u = solve(ls,op)
  # π1u = interpolate(u_vector,T)
  div_π1u = divergence(π1u)

  # compute divergence of vector field, interpolate into Hdiv
  div_u = divergence(u_vector)
  # a_l2(p,q) = ∫(p*q)dΩ
  # l_l2(q) = ∫(div_u*q)dΩ
  # op = AffineFEOperator(a_l2,l_l2,D,E)
  # π2_div_u = solve(op)
  π2_div_u = interpolate(div_u,E)

  e = div_π1u-π2_div_u
  e_hdiv = sqrt(sum(∫( e⋅e )dΩ))

  if return_vtk
    writevtk(Ω,dir*"/de_Rham_H1_nref$(lvl)_p$p_fe",
      cellfields=["grad_perp_π0u"=>grap_perp_π0u,"π1u_grad_perp_u"=>π1u_grad_perp_u,
                "error"=>π1u_grad_perp_u-grap_perp_π0u],append=false);
    writevtk(Ω,dir*"/de_Rham_Hdiv_nref$(lvl)_p$p_fe",nsubcells=10,
      cellfields=["div_π1u"=>div_π1u,"π2_div_u"=>π2_div_u,
              "error"=>div_π1u-π2_div_u],append=false);
  end


  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = num_cells(model)
  dxx = 1/sqrt(n)
  output = @strdict e_h1 e_hdiv n dxx p_fe lvl
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("derham_nref$(lvl)_p$p_fe.jld2")), output)

  return e_h1, e_hdiv, false

end
