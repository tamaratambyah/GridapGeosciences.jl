using Gridap
using GridapGeosciences
using Gridap.Geometry, Gridap.Arrays, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Fields
import GridapGeosciences.Geometry: _CCAM_cube_nodes_3d
import GridapGeosciences.Helpers: normal_vec
using FillArrays
using DrWatson
dir = datadir("Extruded_model")
!isdir(dir) && mkdir(dir)

a = π/4 ### for testing, use a = 1.0
npanels = 6

data = [ 1,2,3,4,  9,10,11,12,
         3,4,5,6,  11,12,13,14,
         2,7,4,6,  10,15,12,14,
         8,5,7,6,  16,13,15,14,
         1,8,2,7,  9,16,10,15,
         1,3,8,5,  9,11,16,13  ]
ptr = generate_ptr(3,npanels)
cell_node_ids = Table(data,ptr)

polytopes = fill(HEX,1)
cell_type = fill(1,npanels)
reffes = LagrangianRefFE(Float64,HEX,1)
cell_reffes=[reffes]

alpha_beta_gamma = [
  Point(-a, -a, 0.0)
  Point(a,  -a, 0.0)
  Point(-a,  a, 0.0)
  Point(a,   a, 0.0)
  Point(-a, -a, 1.0)
  Point(a,  -a, 1.0)
  Point(-a,  a, 1.0)
  Point(a,  a, 1.0)
]
cell_vertices_alpha_beta_gamma = fill(alpha_beta_gamma,npanels)
scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(HEX,Gridap.ReferenceFEs.lagrangian,Float64,1)
cell_shape_funs =
   FillArrays.Fill( Gridap.ReferenceFEs.get_shapefuns(scalar_reffe), length(cell_vertices_alpha_beta_gamma) )
cmap = lazy_map(linear_combination,cell_vertices_alpha_beta_gamma,cell_shape_funs)

cube_surface_nodes = _CCAM_cube_nodes_3d(a)
nodes =  [cube_surface_nodes
          cube_surface_nodes]

topo = UnstructuredGridTopology(nodes,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
labels = FaceLabeling(topo)
grid = Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented(),nothing,cmap)
extruded_panel_model = UnstructuredDiscreteModel(grid,topo,labels)
# extruded_panel_model_ref = Adaptivity.refine(extruded_panel_model)

writevtk(Triangulation(extruded_panel_model),dir*"/extruded_cube",append=false)

############ refinement
# model = extruded_panel_model
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(1,1)))
method = Adaptivity.RedGreenRefinement()
cells_to_refine = nothing

ctopo = get_grid_topology(model)
coarse_labels = get_face_labeling(model)
# Create new model
rrules, faces_list = Adaptivity.setup_edge_based_rrules(method,ctopo,cells_to_refine)
topo   = Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)
reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))

cell_vertices_alpha_beta_gamma = fill(alpha_beta_gamma,num_cells(topo))
scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(HEX,Gridap.ReferenceFEs.lagrangian,Float64,1)
cell_shape_funs =
   FillArrays.Fill( Gridap.ReferenceFEs.get_shapefuns(scalar_reffe), length(cell_vertices_alpha_beta_gamma) )
cmap = lazy_map(linear_combination,cell_vertices_alpha_beta_gamma,cell_shape_funs)
Dc = num_cell_dims(model)

p = ctopo
D = num_cell_dims(p)
@check length(faces_list) == D+1

coords_old = get_cell_coordinates(get_grid(model))
coords_new = Vector{eltype(coords_old)}(undef,num_cells(topo))

if !isempty(faces_list[1])
  for node in faces_list[1]
    println(node)
    # coords_new[n] = coords_old[node]
    # n += 1
  end
end

d = 1
dfaces = faces_list[2]
d2n_map = get_faces(p,d,0)
cache = array_cache(d2n_map)
for face in dfaces
  face_nodes = getindex!(cache,d2n_map,face)
  println(face_nodes)
  # coords_new[n] = sum(coords_old[face_nodes])/length(face_nodes)
  # n += 1
end

  nN_new     = sum(length,faces_list)
  coords_old = get_vertex_coordinates(p)
  coords_new = Vector{eltype(coords_old)}(undef,nN_new)

  n = 1
  # Nodes
  if !isempty(faces_list[1])
    for node in faces_list[1]
      coords_new[n] = coords_old[node]
      n += 1
    end
  end
  # Faces (d > 0)
  for (d,dfaces) in enumerate(faces_list[2:end])
    if !isempty(dfaces)
      d2n_map = get_faces(p,d,0)
      cache = array_cache(d2n_map)
      for face in dfaces
        face_nodes = getindex!(cache,d2n_map,face)
        coords_new[n] = sum(coords_old[face_nodes])/length(face_nodes)
        n += 1
      end
    end
  end

  return coords_new





grid   = UnstructuredGrid(
  get_vertex_coordinates(topo),
  get_faces(topo,Dc,0),reffes,
  get_cell_type(topo),
  OrientationStyle(topo),
  nothing,
  cmap
)
glue = Adaptivity.blocked_refinement_glue(rrules)
labels = Adaptivity.refine_face_labeling(coarse_labels,glue,ctopo,topo)
ref_model = UnstructuredDiscreteModel(grid,topo,labels)
extruded_panel_model_ref =  AdaptedDiscreteModel(ref_model,model,glue)

writevtk(Triangulation(extruded_panel_model_ref),dir*"/extruded_cube_ref",append=false)

############### map to extruded panel
for p in collect(1:6)
  mask = panel_ids .== p
  writevtk(Triangulation(extruded_panel_model,mask),dir*"/aligned_panel_$p",append=false)
end

include("forward_map_serial.jl")
cell_geo_map = geo_map_func_3D(get_panel_ids(extruded_panel_model_ref))
writevtk(Triangulation(extruded_panel_model_ref),dir*"/extruded_panel_model",append=false,geo_map=cell_geo_map)
