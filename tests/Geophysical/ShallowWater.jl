"""
solve the non-linear shallow water equations in steady form using manufactured solutions
u + q F^⟂ + ∇ᵧ(Φ) = f₁
φ + ∇ᵧ⋅F = f₁
F = φu
Φ = 0.5(u⋅u) + gᵣφ
q = 1/φ( ∇ᵧ^⟂⋅u  + f )
"""

module ShallowWater

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using Test

include("../convergence_tools.jl")
include("Williamson2Test.jl")


function nonlinear_shallow_water_solver(
  panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,dir::String,h::Function,vX::Function,f::Function,η::Function,ls=LUSolver(),
    return_vtk=false,check_geo_balance=false)

  ranks = get_ranks(panel_model)

  lvl = nref(nc(panel_model))
  i_am_main(ranks) &&  println("Refinement level: $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  R = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  H = TrialFESpace(R)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)


  ## initial conditions
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_contra_h = interpolate(u_contra_cf,U)
  # u_proj_h = covarient_basis_cf ⋅ u_contra_h
  u_proj_h = covarient_basis_cf ⋅ u_contra_cf

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  h_h = interpolate(h_cf,P)

  cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  gravity = _g

  # absolute vorticity
  η_cf = panelwise_cellfield(η,Ω_panel,panel_ids)
  η_h = interpolate(η_cf,H)

  # mectrics required in weak forms
  detg_cf = panelwise_cellfield(detg,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)


  #### DIAGNOSTIC VARIABLES
  # mass flux
  biformF(F,v) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ
  # liformF(v) = ∫( h_h*(u_contra_h⋅(metric_cf⋅v))*meas_cf   )dΩ
  liformF(v) = ∫( h_cf*(u_contra_cf⋅(metric_cf⋅v))*meas_cf   )dΩ
  op = AffineFEOperator(biformF,liformF,U,V)
  Fh = solve(ls,op)

  # Bernoulli potential
  biformΦ(Φ,r) = ∫( Φ*r*meas_cf  )dΩ
  # liformΦ(r) = ∫( gravity*h_h*r*meas_cf  )dΩ + ∫( 0.5*( u_contra_h ⋅(metric_cf⋅u_contra_h) )r*meas_cf  )dΩ
  liformΦ(r) = ∫( gravity*h_cf*r*meas_cf  )dΩ + ∫( 0.5*( u_contra_cf ⋅(metric_cf⋅u_contra_cf) )r*meas_cf  )dΩ
  op = AffineFEOperator(biformΦ,liformΦ,P,Q)
  Φh = solve(ls,op)

  # vorticity
  perp_matrix_cf = panelwise_cellfield(perp_matrix,Ω_panel,panel_ids)
  # biformq(q,r) = ∫( q*h_h*r*meas_cf  )dΩ
  biformq(q,r) = ∫( q*h_cf*r*meas_cf  )dΩ
  # liformq(r) = ∫( cor_cf*r*meas_cf  )dΩ + ∫( (perp_matrix_cf⋅u_contra_h)⋅∇(r)  )dΩ
  liformq(r) = ∫( cor_cf*r*meas_cf  )dΩ + ∫( (perp_matrix_cf⋅u_contra_cf)⋅∇(r)  )dΩ
  op = AffineFEOperator(biformq,liformq,H,R)
  qh = solve(ls,op)

  # e_η = l2((η_h - qh*h_h )*meas_cf,dΩ)
  e_η = l2((η_cf - qh*h_cf ),meas_cf,dΩ)

  #### PROGNOSTIC VARIABLES

  # equation for depth:
  # rhs_h = h_h + 1/meas_cf*( Fh⋅grad_meas_cf + meas_cf*(∇⋅Fh)   )
  rhs_h = h_cf + 1/meas_cf*( Fh⋅grad_meas_cf + meas_cf*(∇⋅Fh)   )
  biform_p(p,r) = ∫( (p*r)*meas_cf )dΩ
  liform_p(r) = ( ∫( (rhs_h*r)*meas_cf )dΩ
                - ∫( r*(Fh⋅grad_meas_cf + meas_cf*(∇⋅Fh) )  )dΩ
                )
  op = AffineFEOperator(biform_p,liform_p,P,Q)
  ph = solve(ls,op)
  e_p = l2((h_cf - ph),meas_cf,dΩ) # error in depth


  # equation for velocity
  Fperph = 1/meas_cf*(perp_matrix_cf⋅Fh)
  # rhs_u = u_contra_h + qh*Fperph + (  inv_metric_cf⋅gradient(Φh) )
  rhs_u = u_contra_cf + qh*Fperph + (  inv_metric_cf⋅gradient(Φh) )

  # check geostropohic balance
  geo_balance = qh*Fperph + (  inv_metric_cf⋅gradient(Φh) )
  e_geo_balance = sum(∫( geo_balance )dΩ)
  if check_geo_balance
    @check vector_length(e_geo_balance) < 1e-12
    i_am_main(ranks) && println("Global geostropohic balance error: $e_geo_balance")
  end

  # solve for velocity
  biform_u(u,v) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ
  liform_u(v) = ( ∫( rhs_u⋅(metric_cf⋅v)*meas_cf )dΩ
                  + ∫( Φh*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                  - ∫( qh*( (perp_matrix_cf⋅Fh) ⋅(metric_cf ⋅v))   )dΩ
                    )
  op = AffineFEOperator(biform_u,liform_u,U,V)
  uh = solve(ls,op)

  uh_proj = covarient_basis_cf ⋅ uh

  # e_u = l2( ( uh-u_contra_h  ),meas_cf,dΩ  )
  e_u = l2( (u_proj_h - uh_proj),meas_cf,dΩ) # error in physical velocity u

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = geo_map_func(Ω_panel)
    panel_cfs = [ph, h_cf,ph-h_cf,
                uh_proj, u_proj_h, uh_proj-u_proj_h,
                 ]
    labels = ["ph","p","ep",
              "uh","u","eu",
                ]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  i_am_main(ranks) && println(e_u, "; ", e_p, "; ",  e_η)

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict e_u e_p e_η n dxx p_fe lvl
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("shallow_water_nref$(lvl)_p$p_fe.jld2")), output)

  return e_u,  e_p, e_η

end

function nonlinear_shallow_water_errors(panel_model,p_fe::Int,dir::String,
  h::Function,vX::Function,f::Function,η::Function,b::Function,
 ls=LUSolver(),CFL=0.1,return_vtk=false,check_geo_balance=false)
  nonlinear_shallow_water_solver(panel_model,p_fe,dir,h,vX,f,η,ls,return_vtk,check_geo_balance)
end

################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: ShallowWater")

  n_ref_lvls = 4
  ps = [1,2,3]
  ζs = [0.0]
  ls = LUSolver()
  models  = get_refined_models(n_ref_lvls)

  dir = datadir("ShallowWaterConvergence")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    if octree
      i_am_main(ranks) && println("Octrees")
      models =  get_octree_refined_models(ranks,n_ref_lvls)
    else
      models,  = get_distributed_refined_models(ranks,nprocs,models)
    end
  end

  for (i,ζ) in enumerate(ζs)
    _dir = dir*"/func_z$i"
    (i_am_main(ranks) && !isdir(_dir) ) && mkdir(_dir)

    h = panel_to_cartesian(h₀(ζ))
    vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
    f = panel_to_cartesian(f₀(ζ))
    η = panel_to_cartesian(η₀(ζ))

    i_am_main(ranks) && println("wave_equation_convergence_func_z$i")
    p_convergence_test(ranks,ps,models,nonlinear_shallow_water_solver,_dir,h,vX,f,η,ls,true)
  end

  i_am_main(ranks) && println("--DONE--")

end

################################################################################
#### Convergence test with plots
################################################################################

function nonlinear_shallow_water_convergence_test(ranks::AbstractArray,nprocs::Int,
  ζs=[0.0],n_ref_lvls=4,ps=[1],ls=LUSolver(),CFL=0.1,return_vtk=false,check_geo_balance=false)

  williamson2_convergence_test(ranks,nprocs,nonlinear_shallow_water_errors,ζs,n_ref_lvls,ps,ls,CFL,return_vtk,check_geo_balance)
end


end # module
