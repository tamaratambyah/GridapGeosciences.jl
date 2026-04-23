"""
solve the linearised shallow water equations in steady form using manufactured solutions
u + f (k×u) + ∇ᵧ(φ) = f₁
φ + ∇ᵧ⋅u = f₁
"""

module LinearisedShallowWaterTests

using Gridap
using Gridap.Helpers
using Gridap.Algebra
using GridapGeosciences
using GridapP4est
using Test

# using DrWatson
# using MPI
# using PartitionedArrays

# using GridapPETSc
# function petsc_mumps_setup(ksp)
#   pc       = Ref{GridapPETSc.PETSC.PC}()
#   mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
#   @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
#   @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
#   @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
#   @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
#   @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
#   @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
#   @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 1)
#   @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
#   @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 1)
#   # @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[],  1, 0.00001)
#   @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
# end

include("Williamson_functions.jl")


a_e = 6.37e6 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 2.94e4/g #m
u_0 = 2*π*a_e/T #m/s

L = a_e
_τ = 1/ω

_a = a_e/L
_g = g*_τ^2/L
_ω = ω*_τ
_H_0 = H_0/L
_T = T/_τ
_u0 = u_0/L*_τ


function linear_shallow_water_solver(panel_model,
  p_fe::Int,dir::String,h::Function,vX::Function,f::Function,ls=LUSolver(),return_vtk=false;
  _i_am_main=true)

  Dc = num_cell_dims(panel_model)
  lvl = nref(panel_model)

  _i_am_main && println("nref = $lvl; p_fe = $p_fe; Dc = $Dc")

  degree = 5*(p_fe+1)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_panel,2*degree)

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  # metric information
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω_panel,panel_ids)

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(piola(vX),Ω_panel,panel_ids)
  u_proj_cf = covariant_basis_cf ⋅(1/meas_cf * u_cf  )
  cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)

  p_int = interpolate(h_cf,P)
  u_int = interpolate(u_cf,U)

  ## Here we construct the coriolis term in 2D v. 3D
  ## On the surface, the term is: ∫( ̃f ( ̃k × ̃u  )  )dΩ
  ##
  ## In 2D, we use the rotation matrx
  function vecPerp(u)
    # u   = (u1, u2)
    # u^T = (-u2, u1)
    VectorValue(-u[2],u[1])
  end
  Aperp = [0 -1
          1 0]
  Rperp = TensorValue(Aperp)
  Rperp_cf = CellField(Rperp,Ω_panel)

  ## In 3D, we construct ̃k using the area measure
  _area_meas(p) = x->  forward_jacobian_3D(p,x) ⋅ (inv_metric(p,x) ⋅ VectorValue(1,0,0))
  area_meas(p) = x-> norm(_area_meas(p)(x))
  normal_3D(p) = x-> (1/area_meas(p)(x) )*VectorValue(1,0,0)
  normal_3D_cf = panelwise_cellfield(normal_3D,Ω_panel,panel_ids)

  ## return the appropriate term based on Dimension
  function get_coriolis_term(Dc::Int)
    if Dc == 2
      return ((u,p),(v,q)) -> ∫( ( cor_cf*( (Rperp_cf⋅ u)⋅v))  )dΩ
    elseif Dc == 3
      return ((u,p),(v,q)) -> ∫( cor_cf*( normal_3D_cf ×( metric_cf⋅u*(1/meas_cf)  ) )⋅(metric_cf⋅v)*(1/meas_cf)  )dΩ
    end
  end
  coriolis_term((u,p),(v,q)) = get_coriolis_term(Dc)((u,p),(v,q))

  ## construct bilinear form using coriolis_term
  biform_u((u,p),(v,q)) = ( ∫( (u⋅ (metric_cf⋅v))*(1/meas_cf) )dΩ
                        + coriolis_term((u,p),(v,q))
                         - ∫( p*(∇⋅v) )dΩ
                          )
  biform_p((u,p),(v,q)) = ∫( (p*q)*meas_cf )dΩ + ∫( q*(∇⋅u) )dΩ
  biformX((u,p),(v,q)) = biform_u((u,p),(v,q)) + biform_p((u,p),(v,q))


  # manufacture rhs functions
  function get_liform(Dc::Int)

    # the manufactured solution is exactly the LHS operator
    _liformX((v,q)) = (
      ∫( (u_int⋅ (metric_cf⋅v))*(1/meas_cf) )dΩ
    + ∫( gradient(p_int)⋅v )dΩ # assume regularity to IBP
    + coriolis_term((u_int,p_int),(v,q)) # coriolis term
    + ∫( (p_int*q)*meas_cf )dΩ
    + ∫( q*(∇⋅u_int) )dΩ
    )

    if Dc == 2
      return v -> _liformX(v)
    elseif Dc == 3
      # in 3D, account for the boundary term from IBP
      Γ = BoundaryTriangulation(panel_model;tags=["bottom_boundary","top_boundary"])
      dΓ = Measure(Γ,degree)
      nΓ = get_normal_vector(Γ)
      boundary((v,q)) = ∫( (v⋅nΓ)*p_int )dΓ
      return v -> _liformX(v) - boundary(v)
    end
  end


  op = AffineFEOperator(biformX,get_liform(Dc),X,Y)
  uh,ph = solve(ls,op)

  uh_proj = covariant_basis_cf ⋅ (1/meas_cf*uh)

  _e = u_cf - uh
  e_u =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*(1/meas_cf) )dΩ_error))

  _e = h_cf - ph
  e_p = sqrt(sum(∫( (_e*_e)*meas_cf )dΩ_error))

  if return_vtk
    panel_cfs = [ph, uh_proj, uh_proj-u_proj_cf,ph-h_cf]
    labels = ["p","u_proj","eu","ep"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$(p_fe)_D$Dc",
          cellfields=cellfields,append=false,geo_map=geo_map_func(Ω_panel))
  end

  return e_u, e_p, false

end

################################################################################
#### Launch linearised shallow water equation -- on gadi
################################################################################
# function launch_linearised_shallow_water(ranks,Dc,n_ref,p_fe::Int,dir::String,return_vtk=1)

#   i_am_main(ranks) && println("--START--")
#   i_am_main(ranks) && println("Wave equation: Dc = $Dc")

#   dir_convergence = dir*"/convergence"
#   (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

#   # ensure no MPI task tries to generate the file before the main MPI task has
#   # created the folder
#   PartitionedArrays.barrier(ranks)

#   ζ = 0.0
#   h = panel_to_cartesian(h₀(ζ))
#   vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
#   f = panel_to_cartesian(f₀(ζ))

#   omodel = if Dc == 2
#     ParametricOctreeDistributedDiscreteModel(ranks, radius;
#     num_initial_uniform_refinements=n_ref)
#   elseif Dc == 3
#     Parametric3DOctreeDistributedDiscreteModel(ranks,radius,thickness;
#         num_horizontal_uniform_refinements=n_ref,
#         num_vertical_uniform_refinements=n_ref);
#   end

#   panel_model = omodel.parametric_dmodel


#   GridapPETSc.Init()
#   ls = PETScLinearSolver(petsc_mumps_setup)

#   e_u,e_p, = linear_shallow_water_solver(panel_model,p_fe,dir,h,vX,f,ls,Bool(return_vtk);_i_am_main=i_am_main(ranks))

#   i_am_main(ranks) && println("eu = $e_u, e_p = $e_p")

#   ## convergence output for DrWatson
#   n = nc(panel_model)
#   dxx = dx(panel_model)
#   output = @strdict e_u e_p n dxx p_fe n_ref Dc
#   i_am_main(ranks) && safesave(datadir(dir_convergence, ("linearised_sw_nref$(n_ref)_p$(p_fe)_D$Dc.jld2")), output)


#   GridapPETSc.Finalize()
#   GridapPETSc.gridap_petsc_gc()

#   i_am_main(ranks) && println("--DONE--")

# end



################################################################################
#### Auto convergence test
################################################################################
function main(models::AbstractArray;ps=[2],_i_am_main=true)
  h = panel_to_cartesian(h₀(0.0))
  vX = panel_to_cartesian(tangent_vec(u₀(0.0)))
  f = panel_to_cartesian(f₀(0.0))

  ls = LUSolver()
  dir = @__DIR__
  p_convergence_auto_test(ps,models,linear_shallow_water_solver,dir,h,vX,f,ls;_i_am_main=_i_am_main)
end

# function main(distribute,nprocs;)
#   ranks = distribute(LinearIndices((nprocs,)))

#   n_ref_lvls = 4
# radius = 1
# thickness = 0.19
#   ## Distributed model: 2D
#   models = get_distributed_refined_models(ranks,nprocs,n_ref_lvls,radius)
#   main(models;_i_am_main=i_am_main(ranks))

#   ### P4test model: 2D
#   models = get_octree_refined_models(ranks,n_ref_lvls,radius)
#   main(models;_i_am_main=i_am_main(ranks))

#   ### P4test model: 3D
#   models = get_3D_octree_refined_models(ranks,n_ref_lvls-1,radius,thickness)
#   main(models;ps=[1],_i_am_main=i_am_main(ranks))

# end



end # module
