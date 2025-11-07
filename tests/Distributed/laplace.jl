using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using Gridap.Algebra
using GridapDistributed

using DrWatson

include("missing_overloads.jl")

function laplace_beltrami_solver_3D(panel_model,p_fe::Int,dir::String,f::Function,ls=LUSolver(),return_vtk=false)

  das =  FullyAssembledRows()

  ranks = get_ranks(panel_model)

  i_am_main(ranks) && println("Assembly strategy: $das")

  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

  tags = ["bottom_boundary",  "top_boundary"]

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(das,panel_model)
  dΩ = Measure(Ω_panel,2*p_fe+1)

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1,
                dirichlet_tags=tags)
  U = TrialFESpace(V,f_panel_cf)

  # i_am_main(ranks) && println("Zeromean: ", sum(∫(f_panel_cf*meas_cf)dΩ))
  rhs_cf = - slap_panel_cf

  poisson_biform(u,v) =  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
  poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
  assem = SparseMatrixAssembler(U,V,das)
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V,assem)

  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(U,x)

  e = l2(f_panel_cf-uh,dΩ)

  if return_vtk
    _Ω_panel = Triangulation(panel_model)
    ## call geo_map_func on the panel ids that includes ghost+owned
    cell_geo_map = geo_map_func(get_panel_ids(_Ω_panel))
    panel_cfs = [f_panel_cf,uh,f_panel_cf-uh]
    labels = ["u","uh","eu"]
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
  output = @strdict e n n_h n_v dxx p_fe lvl_h lvl_v
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("laplace_beltrami_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe.jld2")), output)

  return e, false,false
end

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

include("../convergence_tools.jl")
dir = datadir("Laplace3D")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

n_ref_v = 4
n_ref_h = 5

p_fe = 1
ls = LUSolver()
return_vtk = true
fXYZ(XYZ) =  XYZ[1]*XYZ[2]*XYZ[3]
f = panel_to_cartesian(fXYZ)

for v_lvl in n_ref_v:-1:1
  models = get_3D_octree_refined_models(ranks,n_ref_h,v_lvl)
  errors,ns,dxs,slopes = h_convergence_test(models,laplace_beltrami_solver_3D,p_fe,dir,f,ls,return_vtk)
end



# ### analysis
using DrWatson
using DataFrames
include("../convergence_tools.jl")
dir = datadir("Laplace3D/convergence")
df = collect_results(dir)

nref_v = unique(df.lvl_v)
plot()

for lvl_v in nref_v
  errors = df[(df.lvl_v .== lvl_v ),:e]
  dxs = df[(df.lvl_v .== lvl_v ),:dxx]
  ns = df[(df.lvl_v .== lvl_v ),:n]

  slope = convergence_rate(dxs,errors)
  plot_convergence(errors,ns,dxs,slope;leginf=["u"],colors=[palette(:tab10)[lvl_v],palette(:tab10)[lvl_v] ] )
end
plot!(show=true)
savefig(dir*"/convergence_laplace_3D")
