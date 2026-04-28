using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson

## should return ghost+owned panel ids
MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

dir = datadir("Omodel_2D_refinement")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

# level 0
radius = 1.0
omodel = ParametricOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=3)
panel_model = omodel.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
panel_ids = get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk_with_cell_geomap(cell_geo_map,Ω_panel,dir*"/model0",append=false)

function _adapt_model(ranks,model::ParametricOctreeDistributedDiscreteModel)
  cell_partition=get_cell_gids(model.octree_dmodel)
  panel_model = model.parametric_dmodel
  panel_ids = get_panel_ids(panel_model)

  ref_flags=map(ranks,partition(cell_partition),panel_ids,local_views(panel_model)) do rank,indices,pid,lmodel
    flags=zeros(Cint,length(indices))
    # flags[pid .== 1] .= refine_flag

    cmap = get_cell_map(get_grid(lmodel))
    ref_points = get_cell_ref_coordinates(lmodel)
    coords = lazy_map(evaluate,cmap,ref_points)
    cell_geo_map = lazy_map(p -> ForwardMap(p), pid)
    fi = lazy_map(p->Cartesian2SphericalMap(),pid)
    latlon_cell_geo_map = lazy_map(∘, fi, cell_geo_map)
    θϕ = lazy_map(evaluate,latlon_cell_geo_map,coords)

    for (i,_θϕ) in enumerate(θϕ)
      θ = map(x->x[1],_θϕ)
      ϕ = map(x->x[2],_θϕ)
      if ( any(θ .> 6 ) || any(θ .< 0.2 ) ) && any( abs.(ϕ) .< 0.2 ) && pid[i] == 1
        flags[i] = refine_flag
      end
    end

    return flags
  end
  ref_model, adaptivity_glue = Gridap.Adaptivity.adapt(model,ref_flags)
  return ref_model, adaptivity_glue
end




# level 1
omodel1, = _adapt_model(ranks,omodel)
panel_model = omodel1.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk_with_cell_geomap(cell_geo_map,Ω_panel,dir*"/model1",append=false)

using Test
include("../Operators/darcy_mass_conservation.jl")
include("../Operators/analytic_funcs.jl")


func = analytic_funcs[:XYZ]
scalar_field = true
return_vtk = true
p_fe = 1
mass_conservation(panel_model,p_fe,dir,func,scalar_field,return_vtk)



  @test all(all.(map(x->x.< 1e-12,errors)))

  output = @strdict errors ns dxs slopes ps
  safesave(datadir(_dir, ("convergence.jld2")), output)
  plot_convergence_from_saved(_dir,"convergence",["p"])




function _adapt_model_1(ranks,model::ParametricOctreeDistributedDiscreteModel)
  cell_partition=get_cell_gids(model.octree_dmodel)
  panel_model = model.parametric_dmodel
  panel_ids = get_panel_ids(panel_model)

  ref_flags=map(ranks,partition(cell_partition),panel_ids,local_views(panel_model)) do rank,indices,pid,lmodel
    flags=zeros(Cint,length(indices))
    # flags[pid .== 1] .= refine_flag

    cmap = get_cell_map(get_grid(lmodel))
    ref_points = get_cell_ref_coordinates(lmodel)
    coords = lazy_map(evaluate,cmap,ref_points)
    cell_geo_map = lazy_map(p -> ForwardMap(p), pid)
    fi = lazy_map(p->Cartesian2SphericalMap(),pid)
    latlon_cell_geo_map = lazy_map(∘, fi, cell_geo_map)
    θϕ = lazy_map(evaluate,latlon_cell_geo_map,coords)

    for (i,_θϕ) in enumerate(θϕ)
      θ = map(x->x[1],_θϕ)
      ϕ = map(x->x[2],_θϕ)
      if  (any(θ .== 0.0) || any(θ .== 2*π))  && any( abs.(ϕ) .< 0.1 ) && pid[i] == 1
        flags[i] = refine_flag
      end
    end

    return flags
  end
  ref_model, adaptivity_glue = Gridap.Adaptivity.adapt(model,ref_flags)
  return ref_model, adaptivity_glue
end

omodel2, = _adapt_model_1(ranks,omodel1)
panel_model = omodel2.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk_with_cell_geomap(cell_geo_map,Ω_panel,dir*"/model2",append=false)
