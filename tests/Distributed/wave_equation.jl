using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using PartitionedArrays
using MPI
using GridapP4est

using GridapGeosciences
using GridapGeosciences.Distributed

using Test


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

include("missing_overloads.jl")
include("../convergence_tools.jl")
include("williamson_funcs_3D.jl")


function wave_solver_3D(panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
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
  # Ω_panel = Triangulation(with_ghost,panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  tags = ["bottom_boundary",  "top_boundary"]

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  b_cf = CellField(VectorValue(0.0,0.0,0.0),Ω_panel)
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv, dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))
  # U = TrialFESpace(V,b_cf)

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

  ####### solve u
  assem = SparseMatrixAssembler(U,V,das)

  biformU(u,v) =∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ
  liformU(v) = ( ∫(  rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ
                + ∫( h_h*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ )
  op = AffineFEOperator(biformU,liformU,U,V,assem)

  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(U,x)

  uh_proj = covarient_basis_cf ⋅ uh

  e_u = l2( (u_proj_cf - uh_proj)*meas_cf,dΩ) # error in physical velocity u

  ###### solve h
  biformP(p,q) = ∫( (p*q)*meas_cf )dΩ
  liformP(q) = ( ∫( (rhs_scalar*q)*meas_cf )dΩ
                - ∫( q*(u_contra_h⋅grad_meas_cf + meas_cf*(∇⋅u_contra_h) )  )dΩ )
  op = AffineFEOperator(biformP,liformP,P,Q)

  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  ph = FEFunction(P,x)

  e_p = l2((h_cf - ph)*meas_cf,dΩ) # error in depth

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
  n_v = nc_vertical(panel_model)
  dxx = dx(n_h)
  output = @strdict e_u e_p n n_h n_v dxx p_fe lvl_h lvl_v
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("wave_equation_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe.jld2")), output)

  return e_u,e_p,false
end


dir = datadir("Wave3D")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

n_ref_v = 3
n_ref_h = 4

p_fe = 1
ls = LUSolver()
return_vtk = true
ζ = 0.0
h = panel_to_cartesian(h₀(ζ))
vX = panel_to_cartesian(tangent_vec(u₀(ζ)))

for v_lvl in n_ref_v:-1:1
  models = get_3D_octree_refined_models(ranks,n_ref_h,v_lvl)
  errors,ns,dxs,slopes = h_convergence_test(models,wave_solver_3D,p_fe,dir,h,vX,ls,return_vtk)
end




### analysis
using DrWatson
using DataFrames
# using GridapGeosciences
include("../convergence_tools.jl")
dir = datadir("Wave3D/convergence")
# dir = datadir("Wave3D_noghost_proc2/convergence")
df = collect_results(dir)

nref_v = unique(df.lvl_v)

plot()
for lvl_v in nref_v
  e_u = df[(df.lvl_v .== lvl_v ),:e_u]
  e_p = df[(df.lvl_v .== lvl_v ),:e_p]

  dxs = df[(df.lvl_v .== lvl_v ),:dxx]
  ns = df[(df.lvl_v .== lvl_v ),:n]

  slope_u = convergence_rate(dxs,e_u)
  slope_p = convergence_rate(dxs,e_p)
  errors = [e_u;e_p]

  plot_convergence(errors,ns,dxs,slope_u;leginf=["u", "p"],colors=[palette(:tab10)[lvl_v],palette(:tab10)[lvl_v]] )

end
plot!(show=true)
savefig(dir*"/wave_equation_3D")
