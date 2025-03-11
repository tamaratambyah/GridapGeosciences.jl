using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools
include("cube_cell_maps.jl")

struct ManifoldGrid{Dc,Dp,Dp_topo,Tp,O,Tn,A} <: Grid{Dc,Dp} where Tn
  cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn}
  phys_grid::UnstructuredGrid{Dc,Dp,Tp,O,Tn} ## physical domain (for FEM)
  phys_cell_map::A ## map to integration space
  # manifold_cell_map::A ### map to manifold - not the integration space
  manifold_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn} ## ambient space (manifold)
end

function CSGrid(cube_model::DiscreteModel,nodes_function::Function)
  cube_grid = get_grid(cube_model)
  panel_ids = get_panel_ids(cube_model)
  CSGrid(cube_grid,panel_ids,nodes_function)
end

function CSGrid(cube_grid::UnstructuredGrid{Dc,Dp_topo,Tp,O,Tn},
  panel_ids::AbstractArray{<:Int64},nodes_function) where {Dc,Dp_topo,Tp,O,Tn}

  cell_coords = get_cell_coordinates(cube_grid)
  cell_node_ids = get_cell_node_ids(cube_grid)
  cmaps = get_cell_map(cube_grid)

  # phys_coords,  manifold_coords, coords_panel1_2D = nodes_function(panel_ids)

  # coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
  # coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
  # coord_panel1_latlon = lazy_map(GnomonicMap(),coords_panel1_2D)
  # coords_panel1_3D = lazy_map(Sigma(), coord_panel1_latlon)
  # manifold_coords = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)
  # phys_coords = lazy_map(Sigma(), manifold_coords)

  coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
  coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
  central_angles = lazy_map(CentralAngleMap(),coords_panel1_2D)
  coord_panel1_latlon = lazy_map(OtherGnomonicMap(), central_angles)
  coords_panel1_3D = lazy_map(Sigma(), coord_panel1_latlon)
  manifold_coords = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)
  phys_coords = lazy_map(Sigma(), manifold_coords)


  ## evaluate the physical nodes
  T = eltype(eltype(phys_coords))
  phys_nodes = similar(phys_coords, T, num_nodes(cube_grid))
  get_nodes_from_coords!(phys_nodes,cell_node_ids,phys_coords)

  ## define grid of physical space (for FEM)
  phys_grid = Gridap.Geometry.UnstructuredGrid(phys_nodes,cell_node_ids,
      get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid))
  Dp = num_point_dims(phys_grid)

  ## cell map to integration
  phys_cell_map = lazy_map(CellPanelMaps(Rp1,R1p,Bump), panel_ids, cmaps)
  A = typeof(phys_cell_map)


  ## evaluate the manifold nodes
  T = eltype(eltype(manifold_coords))
  manifold_nodes = similar(manifold_coords, T, num_nodes(cube_grid))
  get_nodes_from_coords!(manifold_nodes,cell_node_ids,manifold_coords)

  ## define grid of ambient space (manifold)
  manifold_grid = Gridap.Geometry.UnstructuredGrid(manifold_nodes,cell_node_ids,
      get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid))



  ManifoldGrid{Dc,Dp,Dp_topo,Tp,O,Tn,A}(cube_grid,phys_grid,phys_cell_map,manifold_grid)

end


"""
grid API
"""
function Gridap.Geometry.get_node_coordinates(grid::ManifoldGrid)
  Gridap.Geometry.get_node_coordinates(grid.phys_grid)
end

function Gridap.Geometry.get_cell_node_ids(grid::ManifoldGrid)
  Gridap.Geometry.get_cell_node_ids(grid.phys_grid)
end

function Gridap.Geometry.get_reffes(grid::ManifoldGrid)
  Gridap.Geometry.get_reffes(grid.phys_grid)
end

function Gridap.Geometry.get_cell_type(grid::ManifoldGrid)
  Gridap.Geometry.get_cell_type(grid.phys_grid)
end

function Gridap.Geometry.get_cell_map(grid::ManifoldGrid)
  grid.phys_cell_map
end

function get_manifold_cell_map(grid::ManifoldGrid)
  grid.manifold_cell_map
end

function get_manifold_grid(grid::ManifoldGrid)
  grid.manifold_grid
end

function get_cube_grid(grid::ManifoldGrid)
  grid.cube_grid
end

Gridap.Geometry.OrientationStyle(grid::ManifoldGrid) = Gridap.Geometry.NonOriented()


####
a = 1.0
r = a/sqrt(3)
Sigma() = SigmaMap(r)
PanelMap() = PanelRotationMap(rp1)
InvPanelMap() = PanelRotationMap(r1p)
BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)

cubed_sphere_nodes(x) = 1
phys_grid = CSGrid(cube_model_3D,cubed_sphere_nodes)
phys_cmap = get_cell_map(phys_grid)

test_cell_maps(phys_cmap,ref_cell_coords,cell_coords)


ref_model = Gridap.Adaptivity.refine(cube_model_3D)

ref_grid = CSGrid(ref_model,cubed_sphere_nodes)
test_cell_maps(get_cell_map(ref_grid),get_cell_ref_coordinates(ref_model),get_cell_coordinates(ref_model))

writevtk(get_manifold_grid(ref_grid),dir*"/ref_grid",append=false)

ref_ref_model = Gridap.Adaptivity.refine(ref_model)
ref_ref_grid = CSGrid(ref_ref_model,cubed_sphere_nodes)
test_cell_maps(get_cell_map(ref_ref_grid),get_cell_ref_coordinates(ref_ref_model),get_cell_coordinates(ref_ref_model))

writevtk(get_manifold_grid(ref_ref_grid),dir*"/og_ref_ref_grid",append=false)




model = ref_ref_model
cube_grid = get_grid(model)
cell_coords = get_cell_coordinates(model)
panel_ids = get_panel_ids(model)

coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
og_coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
# og_central_angles = lazy_map(CentralAngleMap(),og_coords_panel1_2D)./1

Ncells_per_panel = Int(length(panel_ids)/6)

panel1_domain = ((-1,1,-1,1)) # atan(-1), atan(1)
partition = (Int(sqrt(Ncells_per_panel)),Int(sqrt(Ncells_per_panel)))
panel1_grid = UnstructuredGrid(CartesianGrid(panel1_domain,partition))

central_angles = get_cell_coordinates(panel1_grid)./1
all_central_angles = repeat(central_angles,6)
# coords_panel1_2D = lazy_map(InverseCentralAngleMap(), central_angles)  # return local cartesian coords (non-uniform)

p1 = findall(x->x==1,panel_ids)
og_coords_panel1_2D[p1] .== central_angles
sum( og_coords_panel1_2D.==all_central_angles )

n = Int(sqrt(length(p1)))
glue = ref_ref_model.glue
n2o_faces_map = glue.n2o_faces_map[end][p1]
n2o_cell_to_child = glue.n2o_cell_to_child_id[p1]

n_old_cells = 2
n_new_cells = 4

[1,2,n_new_cells+1,n_new_cells+2] .+ 2

[2*n_new_cells+1,2*n_new_cells+2,3*n_new_cells+1,3*n_new_cells+2]

old_cell = 1
_old_cell = 1:2
xshift = mod( old_cell + n_old_cells + 1, 4)
[(old_cell-1)*n_new_cells + 1,   (old_cell-1)*n_new_cells + 2,
  old_cell*n_new_cells+1,   old_cell*n_new_cells+2] .+ xshift




# function cubed_sphere_nodes(panel_ids)
  Ncells_per_panel = Int(length(panel_ids)/6)

  panel1_domain = (-π/4,π/4,-π/4,π/4) # atan(-1), atan(1)
  partition = (Int(sqrt(Ncells_per_panel)),Int(sqrt(Ncells_per_panel)))
  panel1_grid = UnstructuredGrid(CartesianGrid(panel1_domain,partition))

  central_angles = get_cell_coordinates(panel1_grid)

  # coords_panel1_2D = lazy_map(InverseCentralAngleMap(), central_angles)  # return local cartesian coords (non-uniform)

  # return lat lons
  # panel1_latlon = lazy_map(OtherGnomonicMap(), central_angles)

  # all_panel1_latlon =  repeat(panel1_latlon,6)
  # panel1_sphere = lazy_map(Sigma(),all_panel1_latlon)
  # panelp_sphere = lazy_map(InvPanelMap(), panel1_sphere, panel_ids)
  # panelp_latlon = lazy_map(Sigma(), panelp_sphere)



  using Plots
  plot_coords(panelp_latlon,panel_ids,"test","lon","lat")

  return panelp_latlon,  panelp_sphere, coords_panel1_2D
# end



  markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
  _colors = palette(:tab10)
  # p1 = plot(title = "Cells")
  # p2 = plot(title = "Cell points")
  p3 = plot(title = "test")

  ids = [1,2,4,3,1] # for plotting Gridap ordered nodes

  cache = array_cache(panelp_sphere)
  for i in eachindex(panelp_sphere)
    println(i)
    panel = panel_ids[i]
    out = getindex!(cache, panelp_sphere, i)

    lon = map(x->x[1],out)
    lat = map(x->x[2],out)
    tt = map(x->x[3],out)

    # plot!(p1,lon[ids],lat[ids],lw=2,c=_colors[panel])
    # scatter!(p2,lon[ids],lat[ids],marker=markers[panel],c=_colors[panel])
    plot!(p3,lon[ids],lat[ids],tt[ids],seriestype=:path,linestyle=:solid,lw=2,
          c=_colors[panel],marker=markers[panel])
  end
  # plot!(p1,legend=false,xlabel="$xlab",ylabel="$ylab")
  # plot!(p2,legend=false,xlabel="$xlab",ylabel="$ylab")
  plot!(p3,show=true)
  # plot!(p3,legend=false,xlabel="$xlab",ylabel="$ylab")

  # savefig(p1,plotsdir()*"/$(simName)_cells")
  # savefig(p2,plotsdir()*"/$(simName)_cells_points")
  # savefig(p3,plotsdir()*"/$(simName)_mesh")

    ## evaluate the manifold nodes
    T = eltype(eltype(panelp_sphere))
    manifold_nodes = similar(panelp_sphere, T, num_nodes(cube_grid))
    get_nodes_from_coords!(manifold_nodes,get_cell_node_ids(cube_grid),panelp_sphere)

    ## define grid of ambient space (manifold)
    manifold_grid = Gridap.Geometry.UnstructuredGrid(manifold_nodes,get_cell_node_ids(cube_grid),
        get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid))


writevtk(manifold_grid,dir*"/ref_ref_grid",append=false)
