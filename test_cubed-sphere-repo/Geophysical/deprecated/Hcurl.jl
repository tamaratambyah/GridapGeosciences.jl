using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc

dir = datadir("SW_3D")
!isdir(dir) && mkdir(dir)

include("../convergence_tools.jl")
include("Williamson2Test.jl")
include(srcdir("Helpers/overloads.jl"))
# include("CurlConformingFESpacesFixes.jl")

## pullback 3D vector to 3D chart
inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = αβ -> inv_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


#################### sphere
n_ref_lvls = 0

o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=n_ref_lvls,
        num_vertical_uniform_refinements=0)
panel_model = o3model.parametric_dmodel

tags = ["top_boundary", "bottom_boundary"]
p_fe = 1

panel_ids = get_panel_ids(panel_model)
Ω = Triangulation(panel_model)
dΩ = Measure(Ω,6)
dΩ_error = Measure(Ω,8*p_fe)

metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω,panel_ids)

# ## normal vector in the chart
# # f_cf = CellField(VectorValue(1,0.0,0.0),Ω)
function fV(p)
  function f(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    VectorValue(xyz[1],xyz[2],xyz[3])
  end
end

################################################################################
#### Consider a more complicated vector field
################################################################################
function fV(p)
  function f(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    VectorValue(-xyz[2],0.0,0.0)
  end
end

f_cf = panelwise_cellfield(contra_v_3D(fV),Ω,panel_ids)


Rcurl = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
RH1 = TestFESpace(panel_model, ReferenceFE(lagrangian,VectorValue{3,Float64},p_fe);conformity=:H1,dirichlet_tags=tags)
RHdiv = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe);conformity=:Hdiv,dirichlet_tags=tags)
RL2 = TestFESpace(panel_model, ReferenceFE(lagrangian,VectorValue{3,Float64},p_fe);conformity=:L2,dirichlet_tags=tags)

R = RH1
H = TrialFESpace(R,f_cf)
f_h = interpolate(f_cf,H)
d = collect(get_cell_dof_values(f_h.fields.item_ref[]))
_d = map(x->round.(x),d)
_d[1]




Rspaces = [Rcurl, RH1, RHdiv, RL2]

# R = Rcurl
for (R,name) in zip(Rspaces, ["curl", "H1", "Hdiv", "L2"])

  H = TrialFESpace(R,f_cf)
  f_h = interpolate(f_cf,H)

  # a(u,v) = ∫( u⋅(metric_cf⋅v)*meas_cf )dΩ
  # l(v) = ∫( f_cf⋅(metric_cf⋅v)*meas_cf )dΩ
  # op = AffineFEOperator(a,l,H,R)
  # A = get_matrix(op)
  # f_h = solve(LUSolver(),op)

  writevtk(Ω,dir*"/sol_$name.vtu",
  cellfields=["f_param"=>f_h, "f_ambient"=>covariant_basis_cf⋅f_h,
              "f"=>f_cf, "famb"=>covariant_basis_cf⋅f_cf,
              "e"=>f_cf-f_h, "eamb"=>covariant_basis_cf⋅(f_cf-f_h) ],
  append=false,geo_map= latlon_geo_map_func(Ω))
end

# a(u,v) = ∫( u⋅(metric_cf⋅v)*meas_cf )dΩ
# l(v) = ∫( f_cf⋅(metric_cf⋅v)*meas_cf )dΩ
# op = AffineFEOperator(a,l,H,R)
# f_h = solve(LUSolver(),op)

grad = gradient(f_cf)
gradh = gradient(f_h)

# ccurl = curl(metric_cf⋅f_cf)
# ccurlh = curl(metric_cf⋅f_h)

eh = f_cf - f_h
_eh = f_cf - f_h

errs[i] =  sum(∫( (_eh⋅_eh)*meas_cf )*dΩ_error)
println("Error lvl $n_ref_lvls: ", sum(∫( (_eh⋅_eh)*meas_cf )*dΩ_error) )
# sum(∫( (eh⋅eh)*meas_cf )*dΩ_error)
# sum(∫( eh⋅eh )*dΩ_error)/sum(∫( f_cf⋅f_cf )*dΩ_error)

eh_grad = grad-gradh
sum(∫( (eh_grad⊙eh_grad)*meas_cf )*dΩ_error)

# eh_curl = ccurl-ccurlh
# sum(∫( (eh_curl⊙eh_curl)*meas_cf )*dΩ_error)

# latlon_geo_map = geo_map_func(Ω)
latlon_geo_map = latlon_geo_map_func(Ω)
panel_cfs = [f_h, f_cf, eh, _eh,
            gradh, grad, eh_grad   ]
cellfields = map((x,y) -> x=>y, ["f_h", "f", "ef_ambient", "ef", "grad_h", "grad", "egrad" ],panel_cfs)
writevtk(Ω,dir*"/sol.vtu", cellfields=cellfields,append=false,geo_map=latlon_geo_map)
