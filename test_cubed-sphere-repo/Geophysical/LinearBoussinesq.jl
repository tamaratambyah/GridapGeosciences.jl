"""
solve the linear Boussineq equations in 3D in steady form using manufactured solutions
u + f(̂n×u) + ∇ᵧ(φ) - bn̂ = f₁
φ + c² ∇ᵧ⋅u = f₂
b + N² u⋅̂n = f₃
"""

module LinearisedBoussinesqTests

using Gridap
using Gridap.Helpers
using Gridap.Algebra
using GridapDistributed
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


import GridapGeosciences.Helpers: RADIUS, THICKNESS
THICKNESS
RADIUS

a_e = 6.37e6/125 # m
d = 5000 #m
Lz = 20e3 #m
R = a_e # m radius
u_0 = 20 #m/s
Ωr = 7.292e-5 #1/s
c = 343 #m/s speed of sound
N = 0.01 #1/s bouyancy frequency
ztop = 10e3 #m
dΘ = 1 #K

LH = a_e # m
LV = ztop/THICKNESS
τ = 1/Ωr # s

_d = d/LV
_Lz = Lz/LV
_R = R/LH
_u_0 = u_0*τ/LH
_Ωr = Ωr*τ
_c = c*τ/LH
_N = N*τ
_ztop = ztop/LV


function p0(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  sin(ϕ)
  #  = 0.0
end

function b0(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr

  θc = 2*π/3
  ϕc = 0.0

  k = sqrt(x^2 + y^2 + z^2) - _R

  r = _R*acos( sin(ϕc)*sin(ϕ) + cos(ϕc)*cos(ϕ)*cos(θ-θc)    )
  s = _d^2/(_d^2 + r^2)
  b = dΘ*s*sin( 2*π*k/_Lz  )
  b
end

function u0(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr

  # u = _u_0*cos(ϕ) #
  # v = 0.0#
  u = -_u_0*y/_R
  v = _u_0*x/_R

  VectorValue(u,v,0.0)
end

function omega(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  2*_Ωr*sin(ϕ)
end



function linear_boussineseq(panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,
  h::Function,vX::Function,f::Function,b::Function,
  ls=LUSolver(),return_vtk=false;
  _i_am_main=true)

  Dc = num_cell_dims(panel_model)
  lvl = nref(panel_model)

  _i_am_main && println("nref = $lvl; p_fe = $p_fe; Dc = $Dc")

  degree = 4*(p_fe+1)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_panel,2*degree)

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe))
  U = TrialFESpace(V)

  W = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe))
  B = TrialFESpace(W)

  Y = MultiFieldFESpace([V, Q, W])
  X = MultiFieldFESpace([U, P, B])


  # metric information
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω_panel,panel_ids)


  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(piola(vX),Ω_panel,panel_ids)
  u_proj_cf = covariant_basis_cf ⋅(1/meas_cf * u_cf  )
  b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)
  omega_cf = panelwise_cellfield(f,Ω_panel,panel_ids)

  p_int = interpolate(h_cf,P)
  u_int = interpolate(u_cf,U)
  b_int = interpolate(b_cf,B)

  ## In 3D, we construct ̃k using the area measure
  _area_meas(p) = x->  forward_jacobian_3D(p,x) ⋅ (inv_metric(p,x) ⋅ VectorValue(1,0,0))
  area_meas(p) = x-> norm(_area_meas(p)(x))
  normal_3D(p) = x-> (1/area_meas(p)(x) )*VectorValue(1,0,0)
  normal_3D_cf = panelwise_cellfield(normal_3D,Ω_panel,panel_ids)

  coriolis_term((u,p,b),(v,q,r)) = ∫( omega_cf*( normal_3D_cf ×( metric_cf⋅u*(1/meas_cf)  ) )⋅(metric_cf⋅v)*(1/meas_cf)  )dΩ
  bouyancy_term(b,v) = ∫( b*(normal_3D_cf⋅v)  )dΩ # v ∈ Hdiv, b ∈ L2

  biform_u((u,p,b),(v,q,r)) = ( ∫( (u⋅ (metric_cf⋅v))*(1/meas_cf) )dΩ
                              + coriolis_term((u,p,b),(v,q,r))
                              - ∫( p*(∇⋅v) )dΩ
                              - bouyancy_term(b,v)
                              )

  #### Pressure
  biform_p((u,p,b),(v,q,r)) = ∫( (p*q)*meas_cf )dΩ + ∫( _c^2*(q*(∇⋅u)) )dΩ

  #### Bouyancy
  biform_b((u,p,b),(v,q,r)) = ∫( (b*r)*meas_cf )dΩ + _N^2* bouyancy_term(r,u)

  biformX((u,p,b),(v,q,r)) = biform_u((u,p,b),(v,q,r)) + biform_p((u,p,b),(v,q,r)) + biform_b((u,p,b),(v,q,r))

  # Account for the boundary term from IBP in the RHS forcing operator
  Γ = BoundaryTriangulation(panel_model;tags=["bottom_boundary","top_boundary"])
  dΓ = Measure(Γ,degree)
  nΓ = get_normal_vector(Γ)

  liformX((v,q,r)) = (
    ∫( (u_int⋅ (metric_cf⋅v))*(1/meas_cf) )dΩ
  + ∫( gradient(p_int)⋅v )dΩ # assume regularity to IBP
  + coriolis_term((u_int,p_int,b_int),(v,q,r)) # coriolis term
  - bouyancy_term(b_int,v)
  + ∫( (p_int*q)*meas_cf )dΩ + ∫( _c^2*(q*(∇⋅u_int)) )dΩ
  + ∫( (b_int*r)*meas_cf )dΩ + _N^2* bouyancy_term(r,u_int)
  - ∫( (v⋅nΓ)*p_int )dΓ
  )


  #### Multifield problem
  op = AffineFEOperator(biformX,liformX,X,Y)
  A = get_matrix(op)
  b_vec = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b_vec)
  xh = FEFunction(X,x)
  uh,ph,bh = xh


  uh_proj = covariant_basis_cf ⋅ (1/meas_cf*uh)

  _e = u_cf - uh
  e_u =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*(1/meas_cf) )dΩ_error))

  _e = h_cf - ph
  e_p = sqrt(sum(∫( (_e*_e)*meas_cf )dΩ_error))

  _e = b_cf - bh
  e_b = sqrt(sum(∫( (_e*_e)*meas_cf )dΩ_error))

  if return_vtk
    panel_cfs = [h_cf, u_proj_cf, b_cf, ph, uh_proj, bh, h_cf-ph, u_proj_cf-uh_proj , b_cf-bh]
    labels = ["p","u_proj", "b", "ph", "uh_proj", "bh", "ep","eu", "eb"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",
    cellfields=cellfields,append=false,geo_map=geo_map_func(Ω_panel))
  end

  return e_u, e_p, e_b
end

################################################################################
#### Launch linearised Boussineq equation -- on gadi
################################################################################
# function launch_linearised_boussineseq(ranks,Dc,n_ref,p_fe::Int,dir::String,return_vtk=1)

#   i_am_main(ranks) && println("--START--")
#   i_am_main(ranks) && println("Linearised Boussineq: Dc = $Dc")

#   dir_convergence = dir*"/convergence"
#   (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

#   # ensure no MPI task tries to generate the file before the main MPI task has
#   # created the folder
#   PartitionedArrays.barrier(ranks)

#   h = panel_to_cartesian(p0)
#   vX = panel_to_cartesian(tangent_vec(u0))
#   f = panel_to_cartesian(omega)
#   b = panel_to_cartesian(b0)

# radius,thickness = 1.0, 0.19
#   Parametric3DOctreeDistributedDiscreteModel(ranks,radius,thickness;
#         num_horizontal_uniform_refinements=n_ref,
#         num_vertical_uniform_refinements=n_ref);
#   panel_model = omodel.parametric_dmodel


#   GridapPETSc.Init()
#   ls = PETScLinearSolver(petsc_mumps_setup)

#   e_u, e_p, e_b = linear_boussineseq(panel_model,p_fe,dir,h,vX,f,b,ls,Bool(return_vtk);_i_am_main=i_am_main(ranks))

#   i_am_main(ranks) && println("eu = $e_u, e_p = $e_p, e_b = $e_b")

#   ## convergence output for DrWatson
#   n = nc(panel_model)
#   dxx = dx(panel_model)
#   output = @strdict e_u e_p e_b n dxx p_fe n_ref Dc
#   i_am_main(ranks) && safesave(datadir(dir_convergence, ("linearised_sw_nref$(n_ref)_p$(p_fe)_D$Dc.jld2")), output)


#   GridapPETSc.Finalize()
#   GridapPETSc.gridap_petsc_gc()

#   i_am_main(ranks) && println("--DONE--")

# end


################################################################################
#### Auto convergence test
################################################################################
function main(models::AbstractArray;ps=[1],_i_am_main=true)
  h = panel_to_cartesian(p0)
  vX = panel_to_cartesian(tangent_vec(u0))
  f = panel_to_cartesian(omega)
  b = panel_to_cartesian(b0)

  ls = LUSolver()
  dir = @__DIR__
  p_convergence_auto_test(ps,models,linear_boussineseq,dir,h,vX,f,b,ls;_i_am_main=_i_am_main)

end

# function main(distribute,nprocs;)
#   ranks = distribute(LinearIndices((nprocs,)))

#   n_ref_lvls = 3
# radius,thickness = 1,0.19
#   ### P4test model: 3D
#   models = get_3D_octree_refined_models(ranks,n_ref_lvls,radius,thickness)
#   main(models;_i_am_main=i_am_main(ranks))

# end


end #module
