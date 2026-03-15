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

dir = datadir("Hcurl")
!isdir(dir) && mkdir(dir)

include("../../convergence_tools.jl")
include("../../Geophysical/Williamson2Test.jl")
include(srcdir("Helpers/overloads.jl"))
# include("../../Geophysical/CurlConformingFESpacesFixes.jl")

## pullback 3D vector to 3D chart
inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = αβ -> inv_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


#################### sphere
# o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
#         num_horizontal_uniform_refinements=0,
#         num_vertical_uniform_refinements=0)
# _panel_model = o3model.parametric_dmodel


# ## normal vector in the chart
# function fV(p)
#   function f(γαβ)
#     xyz = forward_map_3D(p)(γαβ)
#     VectorValue(xyz[1],xyz[2],xyz[3])
#   end
# end

#### Consider a more complicated vector field
function fV(p)
  function f(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    VectorValue(0.0,xyz[3],xyz[1]^2)
  end
end

function nedelec_convergence(_panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,fV::Function,return_vtk)

  lvl = nref(nc_horizontal(_panel_model))

  ### get serial model for simplificity
  panel_model = _panel_model.models.item_ref[]

  tags = ["top_boundary", "bottom_boundary"]

  panel_ids = get_panel_ids(panel_model)
  Ω = Triangulation(panel_model)
  dΩ = Measure(Ω,6)

  mask = panel_ids .== 6
  Ω_p6 = Triangulation(Ω,mask)
  dΩ_p6 = Measure(Ω_p6,8)

  metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω,panel_ids)

  vec_contra_cf = panelwise_cellfield(contra_v_3D(fV),Ω,panel_ids)
  vec_proj_cf = covarient_basis_cf ⋅ vec_contra_cf

  ### lowest order nedelec
  reffe =  ReferenceFE(nedelec,Float64,p_fe)
  R = TestFESpace(panel_model,reffe;conformity=:Hcurl,dirichlet_tags=tags)
  H = TrialFESpace(R,vec_contra_cf)

  ########### Interpolation
  vec_contra_h = interpolate(vec_contra_cf,H)

  _e = vec_contra_cf - vec_contra_h
  el2_interp =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_p6))

  ############ L2 projection
  a(u,v) = ∫( u⋅( metric_cf⋅v)*meas_cf )dΩ
  l(v) = ∫( vec_contra_cf⋅(metric_cf⋅v)*meas_cf )dΩ
  op = AffineFEOperator(a,l,H,R)
  vec_contra_h = solve(LUSolver(),op)

  vec_proj_h = covarient_basis_cf ⋅ vec_contra_h

  _e = vec_contra_cf - vec_contra_h
  el2_proj =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_p6))

  if return_vtk
    panel_cfs = [vec_proj_cf, vec_proj_h, vec_proj_cf-vec_proj_h,
                  ]
    labels = ["u_proj", "u_projh", "eproj" ]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω,dir*"/ambient_model_nref$(lvl)_p$(p_fe)",cellfields=cellfields,
          append=false,geo_map=geo_map_func(Ω))
  end


  return el2_interp, el2_proj,false
end

n_ref_lvls = 4
ps = [0]

dir = datadir("InterpolationConvergence")
!isdir(dir) && mkdir(dir)

models  = get_3D_octree_refined_models(ranks,n_ref_lvls)


_dir = dir*"/vector_func_3D_Hcurl"
!isdir(_dir) && mkdir(_dir)
p_convergence_test(ranks,ps,models,nedelec_convergence,_dir,fV,true)
plot_convergence_from_saved(_dir,"convergence",["Interp", "Proj"])
