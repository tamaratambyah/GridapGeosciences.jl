using DrWatson
# add Gridap#cubed-sphere
using GridapP4est
using Gridap
using MPI
using PartitionedArrays

using GridapDistributed
using GridapGeosciences

include("../convergence_tools.jl")


dir = datadir("derham_trig")
!isdir(dir) && mkdir(dir)

n_ref_lvls = 3

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=1)

function initial_panel_refinement(omodel::ParametricOctreeDistributedDiscreteModel)
  cell_partition=get_cell_gids(omodel.octree_dmodel)
  panel_model = omodel.parametric_dmodel
  panel_ids = get_panel_ids(panel_model)

  ref_flags=map(partition(cell_partition),panel_ids,local_views(panel_model)) do indices,pid,lmodel
    flags=zeros(Cint,length(indices))
    flags.=nothing_flag

    # flags[pid .== 1] .= refine_flag
    flags[pid .== 4] .= refine_flag

    return flags
  end
  return ref_flags
end

ref_coarse_flags = initial_panel_refinement(omodel)
fmodel, adaptivity_glue = Gridap.Adaptivity.adapt(omodel,ref_coarse_flags)


function refine_all(omodel::ParametricOctreeDistributedDiscreteModel)
  cell_partition=get_cell_gids(omodel.octree_dmodel)
  ref_coarse_flags=map(partition(cell_partition)) do indices
    flags=zeros(Cint,length(indices))
    flags.=refine_flag
  end
  fmodel, = Gridap.Adaptivity.adapt(omodel,ref_coarse_flags)
  return fmodel
end


models = Vector{ParametricOctreeDistributedDiscreteModel}(undef,n_ref_lvls+1)
models[end] = fmodel
for (i,n) in enumerate(n_ref_lvls:-1:1)
  fmodel = refine_all(fmodel)
  models[n] = fmodel
end

# # plot all models
# for i in 1:length(models)
#   ref = length(models)-i
#   model = models[i].parametric_dmodel
#   Ω = Triangulation(model)
#   cell_geo_map = latlon_geo_map_func(Ω)
#   writevtk(Ω,dir*"/model_$ref",append=false, geo_map=cell_geo_map);
# end

include("../Geophysical/Williamson2Test.jl")
u(xyz) = VectorValue(-xyz[2],xyz[1],0.0)
vX = panel_to_cartesian(tangent_vec(u))
# vX = panel_to_cartesian(tangent_vec(u₀(0.0)))
# h = panel_to_cartesian(h₀(0.0))


ps = [1,2]
ls = LUSolver()
simName ="convergence"

## test convergence
p_convergence_test(ranks,ps,models,Hdiv,dir,vX,ls,false,true)
plot_convergence_from_saved(dir,simName,["Interp","div"])

## test div
p_convergence_test(ranks,ps,models,Hdiv,dir,vX,ls,true,true)
plot_convergence_from_saved(dir,simName,["sum div","m"])


function Hdiv(
  model::ParametricOctreeDistributedDiscreteModel,
  p_fe::Int,dir::String,vX::Function,ls=LUSolver(),conff=false,return_vtk=false)

  panel_model = model.parametric_dmodel

  lvl = Int(floor(log2(num_cells(panel_model))))

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  V = FESpace(panel_model, ReferenceFE(raviart_thomas, Float64, p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  W = FESpace(panel_model, ReferenceFE(lagrangian, Float64, p_fe); conformity=:L2)
  R = TrialFESpace(W)

  covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)

  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)

  Πdiv_u = interpolate(sdiv_cf,R)
  sum(∫( Πdiv_u  )dΩ )

  u_h = interpolate(u_contra_cf,U)
  div_uh =  (u_h⋅grad_meas_cf + meas_cf*divergence(u_h) )
  sum(∫( div_uh  )dΩ )

  conf = abs(sum(∫( div_uh )dΩ ))

  e = covariant_basis_cf⋅u_contra_cf - covariant_basis_cf⋅u_h
  e_u = sqrt(sum(∫( e⋅e )dΩ))

  e_div = div_uh - Πdiv_u
  e_hdiv = sqrt(sum(∫( e_div⋅e_div )dΩ))

  if return_vtk
    cell_geo_map = latlon_geo_map_func(Ω_panel)
    writevtk(Ω_panel,dir*"/Hdiv_nref$(lvl)_p$p_fe",
      cellfields=["uh"=>covariant_basis_cf⋅u_h,"error_u"=>e,
              "div_u"=>Πdiv_u,"div_uh"=>div_uh,
              "error_div"=>Πdiv_u-div_uh],append=false, geo_map=cell_geo_map);
  end

  if conff
    return  conf, false, false
  else
    return e_u, e_hdiv, false
  end
end
