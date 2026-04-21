""" Poisson problem using Laplace-Beltrami operator
u + Δᵧ(u) = f
Need to remove the kernal via zeromean FE space
"""

module LaplaceBeltramiTests

using Gridap
using Gridap.Helpers
using Gridap.Algebra
using GridapGeosciences
using GridapP4est
using Test


function fX(forward_map)
  function _f(αβ)
    x = forward_map(αβ)
    x[1]*x[2]*x[3]
  end
end


function laplace_beltrami_solver(panel_model,
  p_fe::Int,dir::String,f::Function,ls=LUSolver(),return_vtk=false;
  _i_am_main=true)


  Dc = num_cell_dims(panel_model)
  lvl = nref(panel_model)

  _i_am_main && println("nref = $lvl; p_fe = $p_fe; Dc = $Dc")

  degree = 6*(p_fe+1)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_panel,2*degree)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1, constraint=:zeromean)
  U = TrialFESpace(V)

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)

  @check sum(∫(f_panel_cf*meas_cf)dΩ) < 1e-14 "Function must be zero mean to solve with zeromean FE space!"

  rhs_cf = - slap_panel_cf

  poisson_biform(u,v) =  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ

  # manufacture rhs functions
  function get_liform(Dc::Int)

    # the manufactured solution
    poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ

    if Dc == 2
      return v -> poisson_liform(v)
    elseif Dc == 3
      # in 3D, account for the boundary term from IBP
      Γ = BoundaryTriangulation(panel_model;tags=["bottom_boundary","top_boundary"])
      dΓ = Measure(Γ,degree)
      nΓ = get_normal_vector(Γ)
      f_int = interpolate(f_panel_cf,U)
      boundary(v) = ∫( ( (inv_metric_cf⋅gradient(f_int) )⋅nΓ)*v*meas_cf )dΓ
      return v -> poisson_liform(v) + boundary(v)
    end
  end

  op = AffineFEOperator(poisson_biform,get_liform(Dc),U,V)

  # uh = solve(ls,op)

  ## for pvectors, the ghost may not be in the prange of the get_matrix
  ## This causes issues with GridapSolvers Krylov solvers, in the allocation of x
  ## To avoid, allocate x based on the domain of A
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(U,x)

  _e = f_panel_cf - uh
  e = sqrt(sum(∫( (_e*_e)*meas_cf )dΩ_error))

  if return_vtk
    panel_cfs = [f_panel_cf,uh,f_panel_cf-uh]
    labels = ["u","uh","eu"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",
        cellfields=cellfields,append=false,geo_map=geo_map_func(Ω_panel))
  end


  return e, false,false
end


################################################################################
#### Auto convergence test
################################################################################
function main(models::AbstractArray;ps=[2],_i_am_main=true)

  ls = LUSolver()
  dir = @__DIR__
  p_convergence_auto_test(ps,models,laplace_beltrami_solver,dir,fX,ls;_i_am_main=_i_am_main)
end



end # module
