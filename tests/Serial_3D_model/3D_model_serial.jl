using Gridap
using GridapGeosciences
using Gridap.Geometry, Gridap.Arrays, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Fields
import GridapGeosciences.Geometry: _CCAM_cube_nodes_3d
import GridapGeosciences.Helpers: normal_vec

using DrWatson
dir = datadir("Extruded_model")
!isdir(dir) && mkdir(dir)

extrusion(x) = x +  sqrt(3)*normal_vec(x)
_intrusion(x) = x - sqrt(3)*normal_vec(x)

a = π/4
x = a*Point(1.0,-1.0,-1.0)
ex = extrusion(x)

normal_vec(x)
normal_vec(ex)


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

a = π/4 ### for testing, use a = 1.0
cube_surface_nodes = _CCAM_cube_nodes_3d(a)

extrusion(γ::Float64,x) = x + γ*sqrt(3)*normal_vec(x)
_intrusion(γ::Float64,x) = x - γ*sqrt(3)*normal_vec(x)

extrusion_nodes_3D =  [extrusion.(0.0,cube_surface_nodes)
                      extrusion.(1.0,cube_surface_nodes)]

topo = UnstructuredGridTopology(extrusion_nodes_3D,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
labels = FaceLabeling(topo)
grid = Gridap.Geometry.UnstructuredGrid(extrusion_nodes_3D,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())
extruded_cube_model = UnstructuredDiscreteModel(grid,topo,labels)
extruded_cube_model_ref = Adaptivity.refine(extruded_cube_model)

writevtk(Triangulation(extruded_cube_model),dir*"/extruded_cube",append=false)
writevtk(Triangulation(extruded_cube_model_ref),dir*"/extruded_cube_ref",append=false)

############### map to extruded panel
function extruded_cube_to_αβγ(p::Int)
  if p == 1
    A = [0 1 0
        0 0 1
        1 0 0] # γ is X component
  elseif p == 2
    A = [0 1 0
        -1 0 0
        0 0 1] # γ is Z component
  elseif p == 3
    A = [-1 0 0
          0 0 1
          0 1 0]  # γ is Y component
  elseif p == 4
    A = [0 0 1
         0 1 0
         -1 0 0] # γ is X component
  elseif p == 5
    A = [-1 0 0
          0 1 0
         0 0 -1]  # γ is Z component
  elseif p == 6
    A = [0 0 1
        -1 0 0
         0 -1 0]  # γ is Y component
  end

  return TensorValue(A)

end


panel_normal = [
  Point(1.0, 0.0, 0.0) #+X
  Point(0.0, 0.0, 1.0) #+Z
  Point(0.0, 1.0, 0.0) #+Y
  Point(-1.0, 0.0, 0.0) #-X
  Point(0.0, 0.0, -1.0) #-Z
  Point(0.0, -1.0, 0.0) #-Y
]


A_excube2panel = map(p->extruded_cube_to_αβγ(p),collect(1:6))
b_excube2panel = a*VectorValue(0.0, 0.0, -1.0)

model = extruded_cube_model_ref
panel_ids = get_panel_ids(model)
extruded_cube_cmaps = get_cell_map(get_grid(model))

# extrudedcube2panel_map = lazy_map(p-> MyAffineField(A_excube2panel[p],b_excube2panel), panel_ids)
extrudedcube2panel_map = lazy_map(p-> MyAffineField(A_excube2panel[p],b_excube2panel)∘ Mapp(panel_normal[p]), panel_ids)


extruded_panel_cmaps = lazy_map(∘,extrudedcube2panel_map,extruded_cube_cmaps)
ref_pts = get_cell_ref_coordinates(get_grid(model))
extruded_panel_coords = lazy_map(evaluate,extruded_panel_cmaps,ref_pts)
extruded_cube_coords = lazy_map(evaluate,extruded_cube_cmaps,ref_pts)

cell = 5
cube2 = extruded_cube_coords[cell]

X = cube2[4]
n = panel_normal[1]
γ = X⋅n - a
_intrusion(γ::Float64,x) = x - γ*sqrt(3)*normal_vec(x)
y = _intrusion(γ,X)
_X =  _intrusion(γ,X) + γ*n
αβγ = MyAffineField(A_excube2panel[1],b_excube2panel)(_X)
( MyAffineField(A_excube2panel[1],b_excube2panel)∘ Mapp(panel_normal[1]) )(X)



####### extruded panel model
extruded_cube_grid = get_grid(model)
extruded_topo = get_grid_topology(model)

extruded_panel_nodes = get_node_coordinates(extruded_cube_grid) # these are just junk nodes, never used
extruded_panel_grid = Geometry.UnstructuredGrid(extruded_panel_nodes,
    get_cell_node_ids(extruded_cube_grid),get_reffes(extruded_cube_grid),
    get_cell_type(extruded_cube_grid),OrientationStyle(extruded_cube_grid),
    nothing,extruded_panel_cmaps)
extruded_panel_topo = UnstructuredGridTopology(extruded_panel_nodes,get_cell_node_ids(extruded_cube_grid),
      get_cell_type(extruded_topo),get_polytopes(extruded_topo),OrientationStyle(extruded_topo))
extruded_panel_labels = FaceLabeling(extruded_panel_topo)

extruded_panel_model = ParametricDiscreteModel(extruded_panel_grid,extruded_panel_topo,extruded_panel_labels,panel_ids)

for p in collect(1:6)
  mask = panel_ids .== p
  writevtk(Triangulation(extruded_panel_model,mask),dir*"/aligned_panel_$p",append=false)
end

include("forward_map_serial.jl")
cell_geo_map = geo_map_func_3D(get_panel_ids(extruded_panel_model))
writevtk(Triangulation(extruded_panel_model),dir*"/extruded_panel_model",append=false,geo_map=cell_geo_map)
