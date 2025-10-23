using Gridap
using GridapGeosciences
using Gridap.Geometry, Gridap.Arrays, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Fields
import GridapGeosciences.Geometry: _CCAM_cube_nodes_3d
import GridapGeosciences.Helpers: normal_vec

using DrWatson
dir = datadir("Extruded_model")
!isdir(dir) && mkdir(dir)

function generate_ptr(n)
  nvertices = 8
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

npanels = 6
data = [ 1,2,3,4,5,6,7,8,
         3,4,9,10,7,8,11,12,
         2,13,4,10,6,14,8,12,
         15,9,13,10,16,11,14,12,
         1,15,2,13,5,16,6,14,
         1,3,15,9,5,7,16,11  ]
ptr = generate_ptr(npanels)
cell_node_ids = Table(data,ptr)

polytopes = fill(HEX,1)
cell_type = fill(1,npanels)
reffes = LagrangianRefFE(Float64,HEX,1)
cell_reffes=[reffes]

a = π/4 ### for testing, use a = 1.0
cube_surface_nodes = _CCAM_cube_nodes_3d(a)

extrusion(x) = x + normal_vec(x)

extrusion_nodes_3D = [
  cube_surface_nodes[1]             #node 1
  cube_surface_nodes[2]             #node 2
  cube_surface_nodes[3]             #node 3
  cube_surface_nodes[4]             #node 4
  extrusion(cube_surface_nodes[1])  #node 5
  extrusion(cube_surface_nodes[2])  #node 6
  extrusion(cube_surface_nodes[3])  #node 7
  extrusion(cube_surface_nodes[4])  #node 8
  cube_surface_nodes[5]             #node 9
  cube_surface_nodes[6]             #node 10
  extrusion(cube_surface_nodes[5])  #node 11
  extrusion(cube_surface_nodes[6])  #node 12
  cube_surface_nodes[7]             #node 13
  extrusion(cube_surface_nodes[7])  #node 14
  cube_surface_nodes[8]             #node 15
  extrusion(cube_surface_nodes[8])  #node 16
 ]


topo = UnstructuredGridTopology(extrusion_nodes_3D,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
labels = FaceLabeling(topo)
grid = Gridap.Geometry.UnstructuredGrid(extrusion_nodes_3D,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())
extruded_cube_model = UnstructuredDiscreteModel(grid,topo,labels)
extruded_cube_model_ref = Adaptivity.refine(extruded_cube_model)

writevtk(Triangulation(extruded_cube_model),dir*"/extruded_cube",append=false)
writevtk(Triangulation(extruded_cube_model_ref),dir*"/extruded_cube_ref",append=false)

panel_alignment_axis = a .* [ ## b_panel2cube
  Point(1.0, 0.0, 0.0) #+X
  Point(0.0, 0.0, 1.0) #+Z
  Point(0.0, 1.0, 0.0) #+Y
  Point(-1.0, 0.0, 0.0) #-X
  Point(0.0, 0.0, -1.0) #-Z
  Point(0.0, -1.0, 0.0) #-Y
]
panel_extrusion_component = [1,3,2,1,3,2] # component of 3D point that relates to extrusion variable
## projection a onto b, where b is an axis
projection(a::VectorValue{3},b::VectorValue{3}) = (a⋅b)*b

x = extrusion_nodes_3D[1]
ex = extrusion(x)
p = 1
posneg = sign(panel_alignment_axis[p][panel_extrusion_component[p]])
projection(x,panel_alignment_axis[p])

x - normal_vec(x)
ex - normal_vec(x)

projection(ex,panel_alignment_axis[p])

function extrusion_variable(p::Int,x::VectorValue{3})
   proj = projection(x,panel_alignment_axis[p])
  #  posneg = sign(panel_alignment_axis[p][panel_extrusion_component[p]])
  #  rescale = posneg*( proj - panel_alignment_axis[p] ) # rescale to start at 0.0
    rescale = proj
    ex = rescale[panel_extrusion_component[p]] # get component relative to panel
    if abs(x[panel_extrusion_component[p]]) == a
      return 0.0
    else
      return ex
    end
end

extrusion_variable(p,x)
extrusion_variable(p,ex)


function map_points_to_surface(p::Int,x::VectorValue{3})
  γ = extrusion_variable(p,x)
  posneg = sign(panel_alignment_axis[p][panel_extrusion_component[p]])
  if γ == 0.0
    surf_point = x
  else
    surf_point = x - normal_vec(x)
  end
  array = [i for i in surf_point.data]
  array[panel_extrusion_component[p]] = γ
  VectorValue(array)
end

p = 4
x = a*Point(-1,1,1)

γ = extrusion_variable(p,x)
posneg = sign(panel_alignment_axis[p][panel_extrusion_component[p]])
if γ == 0.0
  surf_point = x
else
  surf_point = x - normal_vec(x)
end
array = [i for i in surf_point.data]
array[panel_extrusion_component[p]] = γ
VectorValue(array)


# Ax + b: multiplication by an invertible matrix and translation of a vector
struct Map2Surface{A}  <: Field
  p::A
end


function Gridap.Arrays.return_cache(f::Map2Surface,cellx::AbstractArray{<:VectorValue{3}})
  x = first(cellx)
  T = typeof(x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::Map2Surface,cellx::AbstractArray{<:VectorValue{3}} )
  p = f.p
  y = cache
  map!(x -> map_points_to_surface(p,x), y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::Map2Surface,x::VectorValue{3})
  T = typeof(x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Map2Surface,x::VectorValue{3})
  p = f.p
  y = cache
  y = map_points_to_surface(p,x)
  return y
end


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
         1 0 0] # γ is X component
  elseif p == 5
    A = [-1 0 0
          0 1 0
         0 0 1]  # γ is Z component
  elseif p == 6
    A = [0 0 1
        -1 0 0
         0 1 0]  # γ is Y component
  end

  return TensorValue(A)

end


model = extruded_cube_model

panel_ids = get_panel_ids(model)
extruded_cube_cmaps = get_cell_map(get_grid(model))
A_extruedcube_2_extrudedpanel = map(p->extruded_cube_to_αβγ(p),collect(1:6))


## map all 3D points-> surface
##    * in this step, store the extrusion variable
## map surface points -> alpha,beta
## add back the extrusion variable
extrudedcube2panel_map = lazy_map(p->MatMultField( A_extruedcube_2_extrudedpanel[p] ) ∘ Map2Surface(p), panel_ids)

extruded_panel_cmaps = lazy_map(∘,extrudedcube2panel_map,extruded_cube_cmaps)
ref_pts = get_cell_ref_coordinates(get_grid(model))
extruded_panel_coords = lazy_map(evaluate,extruded_panel_cmaps,ref_pts)


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

include("forward_map_3D.jl")
cell_geo_map = geo_map_func_3D(get_panel_ids(extruded_panel_model))
writevtk(Triangulation(extruded_panel_model),dir*"/extruded_panel_model",append=false,geo_map=cell_geo_map)

extruded_model = lazy_map(∘,cell_geo_map,extruded_panel_cmaps)
lazy_map(evaluate,extruded_model,ref_pts)

p = 6
pts = extruded_panel_coords[p][5]
α,β,γ = pts

#### compute XYZ point on surface of inner sphere using 2D forward_map
αβ = Point(α,β)
XYZ_surf = forward_map(p,αβ)
XYZ_surf + normal_vec(XYZ_surf)
