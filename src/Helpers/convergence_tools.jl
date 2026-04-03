"""
auto convergence tests
"""

function p_convergence_auto_test(ps::Vector{Int},models::AbstractArray,convergence_func,dir::String,fargs...)

  ranks = get_ranks(testitem(models))
  # nref == refinement level of the mesh
  # 2^nref == number of cells per panel
  # 1/(2^nref) -> to get positive slopes

  lvls = map(x->1/(2^nref(x)),models)
  # lvls = map(x->num_cells(x),models)

  for (i,p_fe) in enumerate(ps)
    errs = h_convergence_auto_test(models,convergence_func,p_fe,dir,fargs...)
    slope = convergence_rate(lvls,errs)
    i_am_main(ranks) && println("slope = $slope")
    @test slope >= p_fe
  end

end

function h_convergence_auto_test(models::AbstractArray,f,p_fe::Int,dir::String,fargs...)
  errs = Float64[]

  for model in models
    e, = f(model,p_fe,dir,fargs...)
    push!(errs,e)
  end

  return errs
end



"""
get_refined_models
returns an array of refined serial models where
  models[1] == most refined model
  models[end] == coarsest model
"""
function get_refined_models(n_ref_lvls::Int,coarse_model=false)

  panel_model = coarse_parametric_model()
  panel_models = Vector{DiscreteModel}(undef,n_ref_lvls)

  for n in n_ref_lvls:-1:1
    panel_model = Gridap.Adaptivity.refine(panel_model)
    panel_models[n] = panel_model
  end

  if coarse_model
    push!(panel_models,coarse_parametric_model())
  end

  panel_models
end


function get_distributed_refined_models(ranks,nprocs,n_ref_lvls::Int,coarse_s_model=false)
  s_models  = get_refined_models(n_ref_lvls,coarse_s_model)
  dmodels, dpanel_ids, owned_panel_ids = get_distributed_refined_models(ranks,nprocs,s_models)
  dmodels
end


function get_distributed_refined_models(ranks,nprocs,s_models::Vector{<:DiscreteModel})

  # get refined models in serial
  spanel_ids = map(m->get_panel_ids(m),s_models)
  s_model_coarse = s_models[end]

  # extract the models and glues in arrays
  models = map(m->Adaptivity.get_model(m),s_models[1:end-1])
  glues = map(m->get_adaptivity_glue(m),s_models[1:end-1])
  if typeof(s_model_coarse) <: AdaptedDiscreteModel
    push!(models,Adaptivity.get_model(s_model_coarse))
  else
    push!(models,s_model_coarse)
  end

  # partition the processors
  part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_model_coarse)) for rank in 1:nprocs]

  # store the partition for each model
  cell_to_part = Vector{Vector{Int32}}(undef,length(models))

  # get the coarse partition
  coarse_cell_to_part = zeros(Int32,num_cells(s_model_coarse))
  for (rank, cells) in enumerate(part_to_cells)
    coarse_cell_to_part[cells] .= rank
  end
  cell_to_part[end] = coarse_cell_to_part

  # get the refine partition based on glue
  for level in length(models)-1:-1:1
    n2o_cells = glues[level].n2o_faces_map[3]
    cell_to_part[level] = cell_to_part[level+1][n2o_cells]
  end

  # construct array of distributed models, distributed panel ids (all + owned only)
  dmodels = Vector{DistributedParametricDiscreteModel}(undef,length(models))
  dpanel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))
  owned_panel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))

  # the coarsest model is the last in the list
  coarse_dmodel = DiscreteModel(ranks,models[end],cell_to_part[end])
  dpanel_ids[end],owned_panel_ids[end] = distributed_panel_ids(coarse_dmodel,spanel_ids[end])
  dmodels[end] = DistributedParametricDiscreteModel(coarse_dmodel,dpanel_ids[end])

  # loop backwards through refinement levels
  for level in length(models)-1:-1:1
    child = DiscreteModel(ranks,models[level],cell_to_part[level])
    parent = dmodels[level+1]
    glue = DistributedAdaptivityGlue(glues[level],parent,child)
    dpanel_ids[level],owned_panel_ids[level] = distributed_panel_ids(child,spanel_ids[level])

    dmodels[level] = DistributedParametricDiscreteModel(child, dpanel_ids[level])
  end

  return dmodels, dpanel_ids, owned_panel_ids

end


#### get array of dmodel - octree models
function get_octree_refined_models(ranks,n_ref_lvls::Int,coarse_model=false)

  dmodels = Vector{DistributedParametricDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:1)
    parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n)
    dmodels[i] = parametric_octree_dmodel.parametric_dmodel
  end


  if coarse_model
    parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
    push!(dmodels,parametric_octree_dmodel.parametric_dmodel)
  end

  return dmodels
end

function get_3D_octree_refined_models(ranks,n_ref_lvls::Int)

  dmodels = Vector{DistributedParametricDiscreteModel}(undef,n_ref_lvls)

  for (i,n) in enumerate(n_ref_lvls:-1:1)
    octree3_model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
                        num_horizontal_uniform_refinements=n,
                        num_vertical_uniform_refinements=n);
    dmodels[i] = octree3_model.parametric_dmodel
  end

  return dmodels

end

"""
nref
returns the level of refinement
  * in 3D, returns the horizontal refinement
"""

nref(nc) = Int(log2(sqrt(nc))) ## level of refinement

function nref(model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}})
  nref(nc(model))
end

function nref(model::GridapDistributed.DistributedDiscreteModel{3,3})
  nref(nc_horizontal(model))
end

## nc = num cells per panel
function nc(model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}})
  num_cells(model)/6
end
function nc(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  nc_horizontal(model) + _nc_vertical(model)
end



## nc = num cells per panel in horizontal
function nc_horizontal(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})

  grid = get_grid(model)
  gids = get_cell_gids(model)

  ## find the number of cells that are on the surface.
  ## i.e. with γ = 0.0
  ## make sure to extract only the owned
  f = map(local_views(grid),partition(gids)) do grid, cids
    cmap = get_cell_map(grid)
    pts = get_cell_ref_coordinates(grid)
    f = lazy_map(evaluate,cmap,pts)
    g = lazy_map(FindSurfaceCells(),f)

    owned_cells = own_to_local(cids)
    sum(g[owned_cells])
  end

  nsurface = sum(f)
  ncells_per_panel = Int(nsurface/6)
  return ncells_per_panel
end

# return square here so vertical is 'like' horitzontal
function nc_vertical(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  n = _nc_vertical(model)
  return Int(n^2)
end

# the actual number of cells in vertical per panel
function _nc_vertical(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  ncells_per_panel = nc_horizontal(model)
  n = num_cells(model)/6
  _n =  n /ncells_per_panel
  return Int(_n)
end

using Gridap.Arrays
using Gridap.Geometry

"""
find the cells at γ = 0.0
"""
struct FindSurfaceCells  <: Map
end

function Gridap.Arrays.return_cache(f::FindSurfaceCells,cellx::AbstractArray{<:VectorValue{3}})
  y = similar(cellx,true)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::FindSurfaceCells,cellx::AbstractArray{<:VectorValue{3}} )
  y = cache
  y = !isempty(findall(x->x[1]==0.0,cellx))
  return y
end

function Gridap.Arrays.return_cache(f::FindSurfaceCells,x::VectorValue{3})
  y = true
  return y
end

function Gridap.Arrays.evaluate!(cache,f::FindSurfaceCells,x::VectorValue{3} )
  y = cache
  y = !isempty(findall(x[1]==0.0))
  return y
end


## element size
function dx(model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}})
  tmp =  4*π*RADIUS^2/num_cells(model)
  sqrt(tmp)
end

function dx(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  horizontal = dx_horizontal(model)
  vertical = dx_vertical(model)
  horizontal*vertical
end

function dx_horizontal(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  horizontal = 4*π*RADIUS^2/(nc_horizontal(model)*6)
  sqrt(horizontal) ## quads so have to sqrt
end

function dx_vertical(model::GridapDistributed.GenericDistributedDiscreteModel{3,3})
  vertical = THICKNESS/_nc_vertical(model)
  vertical ### single layer, so no sqrt
end

"""
get_ranks
"""
function get_ranks(model::DiscreteModel)
  return [true]
end

function get_ranks(dmodel::GridapDistributed.DistributedDiscreteModel)
get_parts(dmodel)
end



function convergence_rate(dxs,errors)
  x = log10.(dxs)
  y = log10.(abs.(errors))
  linreg = hcat(fill!(similar(x), 1), x) \ y
  linreg[2]
end
