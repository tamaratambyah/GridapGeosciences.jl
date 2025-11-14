"""
solve the linearised wave equation in steady form using manufactured solutions
u + ∇ᵧ(φ) = f₁
φ + ∇ᵧ⋅u = f₁
"""

module WaveEquation

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

function wave_solver(
  panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,dir::String,h::Function,vX::Function,ls=LUSolver(),return_vtk=false)

  # das =  FullyAssembledRows()
  das =  SubAssembledRows()

  ranks = get_ranks(panel_model)

  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(das,panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  ### See https://github.com/tamaratambyah/GridapGeosciences.jl/issues/5
  ### use triangulation for FE space
  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)

  sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  sgrad_cf = panelwise_cellfield(sgrad(h),Ω_panel,panel_ids)
  pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)

  # manufacture rhs functions
  rhs_scalar = h_cf + sdiv_cf
  rhs_vector = u_proj_cf + sgrad_cf
  rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component


  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)

  biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  biform2((u,p),(v,q)) = ∫( (p*q)*meas_cf )dΩ + ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

  biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
  liformX((v,q)) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( (rhs_scalar*q)*meas_cf )dΩ

  assem = SparseMatrixAssembler(X,Y,das)
  op = AffineFEOperator(biformX,liformX,X,Y,assem)
  uh,ph = solve(ls,op)

  uh_proj = covarient_basis_cf ⋅ uh

  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,6*p_fe+1)
  e_u = l2( (u_proj_cf - uh_proj),meas_cf,dΩ_error) # error in physical velocity u
  e_p = l2((h_cf - ph),meas_cf,dΩ_error) # error in depth

  if return_vtk
    lvl = nref(nc(panel_model))

    _Ω_panel = Triangulation(panel_model)
    cell_geo_map = geo_map_func(_Ω_panel)

    if das == FullyAssembledRows()
      cell_geo_map = geo_map_func(get_panel_ids(_Ω_panel))
    end

    panel_cfs = [h_cf, u_proj_cf, ph, uh_proj, h_cf-ph, u_proj_cf-uh_proj ]
    labels = ["p","u_proj", "ph", "uh_proj", "ep","eu"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)

  end

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict e_u e_p n dxx p_fe lvl
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("wave_equation_nref$(lvl)_p$p_fe.jld2")), output)

  return e_u,e_p,false
end

### 3D model
function wave_solver(panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,h::Function,vX::Function,ls=LUSolver(),return_vtk=false)

  # das =  FullyAssembledRows()
  das = SubAssembledRows()

  ranks = get_ranks(panel_model)

  i_am_main(ranks) && println("Assembly strategy: $das")

  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

  panel_ids = get_panel_ids(panel_model)

  ### Use no ghost -- with_ghost causes issues in velocity
  Ω_panel = Triangulation(das,panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))
  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,6*p_fe+1)

  tags = ["bottom_boundary",  "top_boundary"]

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  b_cf = CellField(VectorValue(0.0,0.0,0.0),Ω_panel)
  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv, dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))
  # U = TrialFESpace(V,b_cf)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  h_h = interpolate(h_cf,P)

  u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_contra_h = interpolate(u_contra_cf,U)

  sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  sgrad_cf = panelwise_cellfield(sgrad(h),Ω_panel,panel_ids)
  pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)
  u_proj_h = covarient_basis_cf ⋅ u_contra_cf

  # manufacture rhs functions
  rhs_scalar = h_cf + sdiv_cf
  rhs_vector = u_proj_cf + sgrad_cf
  rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component

  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)

  ####### solve as multifield
  biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  biform2((u,p),(v,q)) = ∫( (p*q)*meas_cf )dΩ + ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

  biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
  liformX((v,q)) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( (rhs_scalar*q)*meas_cf )dΩ

  assem = SparseMatrixAssembler(X,Y,das)
  op = AffineFEOperator(biformX,liformX,X,Y,assem)
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  xh = FEFunction(X,x)
  uh,ph = xh

  uh_proj = covarient_basis_cf ⋅ uh

  e_u = l2( (u_proj_cf - uh_proj),meas_cf,dΩ_error) # error in physical velocity u
  e_p = l2((h_cf - ph),meas_cf,dΩ_error) # error in depth

  if return_vtk
    _Ω_panel = Triangulation(panel_model)
    cell_geo_map = geo_map_func(_Ω_panel)

    if das == FullyAssembledRows()
      cell_geo_map = geo_map_func(get_panel_ids(_Ω_panel))
    end

    panel_cfs = [h_cf, u_proj_cf, ph, uh_proj, h_cf-ph, u_proj_cf-uh_proj ]
    labels = ["p","u_proj", "ph", "uh_proj", "ep","eu"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  n_h = nc_horizontal(panel_model)
  n_v = _nc_vertical(panel_model)
  dxx = dx(panel_model)
  output = @strdict e_u e_p n n_h n_v dxx p_fe lvl_h lvl_v
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("wave_equation_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe.jld2")), output)

  return e_u,e_p,false
end


function wave_errors(panel_model,p_fe::Int,dir::String,
  h::Function,vX::Function,f::Function,η::Function,b::Function,ls=LUSolver(),CFL=0.1,return_vtk=false)
  wave_solver(panel_model,p_fe,dir,h,vX,ls,return_vtk)
end


################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false,threedims=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: Wave Equation")

  dir = foldername("WaveConvergence",octree,threedims)
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  n_ref_lvls = 4
  ps = [2]#[1,2]
  ζs = [0.0]
  ls = LUSolver()
  # ls = CGSolver(JacobiLinearSolver();maxiter=3000,rtol=1e-8,verbose=i_am_main(ranks))

  models = get_models(ranks,nprocs,n_ref_lvls;threedims=threedims,octree=octree)

  for (i,ζ) in enumerate(ζs)
    _dir = dir*"/func_z$i"
    (i_am_main(ranks) && !isdir(_dir) ) && mkdir(_dir)

    h = panel_to_cartesian(h₀(ζ))
    vX = panel_to_cartesian(tangent_vec(u₀(ζ)))

    i_am_main(ranks) && println("wave_equation_convergence_func_z$i")
    p_convergence_test(ranks,ps,models,wave_solver,_dir,h,vX,ls,true)
  end

  i_am_main(ranks) && println("--DONE--")
end

################################################################################
#### Convergence test with plots
################################################################################

function wave_convergence_test(ranks::AbstractArray,nprocs::Int,
  ζs=[0.0],n_ref_lvls=4,ps=[1],ls=LUSolver(),CFL=0.1,return_vtk=false)

  williamson2_convergence_test(ranks,nprocs,wave_errors,ζs,n_ref_lvls,ps,ls,CFL,return_vtk)
end

end # module
