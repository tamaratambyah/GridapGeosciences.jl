using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields

struct XTimes2Field{A} <: Field
  A_mat::A
end


function Gridap.Arrays.return_cache(f::XTimes2Field,
  cellx::AbstractArray{<:VectorValue{2}})
  A = f.A_mat
  x = first(cellx)

  T = typeof(A⋅x)
  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end


function Gridap.Arrays.evaluate!(cache,f::XTimes2Field,
  cellx::AbstractArray{<:VectorValue{2}})
  setsize!(cache,size(cellx))
  y = cache.array
  A = f.A_mat
  map!(x -> A⋅x, y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::XTimes2Field,x::VectorValue{2})
  A = f.A_mat
  T = typeof(A⋅x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::XTimes2Field,x::VectorValue{2})
  y = cache
  A = f.A_mat
  y = A.⋅x
  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:XTimes2Field},
  cellx::AbstractArray{<:VectorValue{2}})
  T = typeof(transpose(f.object.A_mat) )
  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:XTimes2Field},
  cellx::AbstractArray{<:VectorValue{2}})
  c, = cache
  setsize!(c,size(cellx))
  y = c.array
  AT = transpose(f.object.A_mat)
  map!(x -> AT, y, cellx)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:XTimes2Field},x::VectorValue{2})
  T = typeof(transpose(f.object.A_mat) )
  y = zero(T)

  return y
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:XTimes2Field},x::VectorValue{2})
  y = cache
  y = transpose(f.object.A_mat)

  return y
end


model1 = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,2)))
model2 = UnstructuredDiscreteModel(CartesianDiscreteModel((0,2,0,2),(2,2)))

Ω1 = Triangulation(model1)
dΩ1 = CellQuadrature(Ω1,2)
quad_pts_1 = get_cell_points(dΩ1)
RT1 = FESpace(model1,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)

Ω2 = Triangulation(model2)
dΩ2 = CellQuadrature(Ω2,2)
quad_pts_2 = get_cell_points(dΩ2)
RT2 = FESpace(model2,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)

u(x) = VectorValue(1.0,1.0)
uh1 = interpolate_everywhere(u,RT1)

A_mat = TensorValue{2,2}(2,0,0,1)

cmaps = map(x-> XTimes2Field(A_mat), 1:num_cells(model1))

mapped_fields = lazy_map(Broadcasting(push_∇),get_data(uh1),cmaps)
cf_mapped = lazy_map(Broadcasting(∘),mapped_fields,cmaps)

cf2 = CellData.GenericCellField(cf_mapped,Ω2,DomainStyle(uh1))
cf2((quad_pts_2))




# do the piola map, then map to model2
Jt = ∇(phi)
Jt_inv = pinvJt(Jt)
det_Jt = meas(Jt)
change = det_Jt*Jt_inv

_RT = FESpace(ambient_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
free_values = zero_free_values(_RT)
s = get_fe_dof_basis(_RT)
cell_vals =  s(cf_mapped)
gather_free_values!(free_values,_RT,cell_vals)












include("src/initialise.jl")


a = π/4
r = 1.0*sqrt(3.0)
global A_bump, B_bump, b_bump = bump_matrics(a)

cube_grid_3D,topo_3D,face_labels_3D,panel_ids = coarse_cube_surface_3D(a)


Dc = num_cell_dims(cube_grid_3D)
Dp_amb = num_point_dims(cube_grid_3D)

cube_cell_coords_3D = get_cell_coordinates(cube_grid_3D)

coords_panel1_3D = lazy_map(Rp1PanelMap3D(), cube_cell_coords_3D, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1_3D)

latlon_panel1 = lazy_map(GnomonicMap(), coords_panel1_2D)
sphere_panel1 = lazy_map(SigmaMap(r),latlon_panel1)
sphere_panelp = lazy_map(R1pPanelMap3D(), sphere_panel1, panel_ids)



ambient_cell_coords = lazy_map(RoundVectorValues(),sphere_panelp)

_cube_grid_3D, = coarse_cube_surface_3D(1.0)
get_cell_coordinates(_cube_grid_3D)


ambient_nodes = get_nodes_from_coords(cube_grid_3D,ambient_cell_coords)

ambient_grid = Gridap.Geometry.UnstructuredGrid(ambient_nodes,get_cell_node_ids(cube_grid_3D),
    get_reffes(cube_grid_3D),get_cell_type(cube_grid_3D),OrientationStyle(cube_grid_3D))

cmaps = get_cell_map(ambient_grid)
evaluate(cmaps[1],Point(1,1))



manifold_grid = ManifoldGrid(cube,cube_grid_3D,panel_ids)

num_point_dims(manifold_grid)
num_cell_dims(manifold_grid)


model = ManifoldDiscreteModel(manifold_grid)
num_point_dims(model)

ref_model = Adaptivity.refine(model)
num_point_dims(ref_model)
num_cell_dims(ref_model)

writevtk(get_ambient_grid(get_grid(ref_model)),dir*"/ref_grid",append=false)

get_panel_ids(ref_model)

ref_ref_model = Adaptivity.refine(ref_model)
num_point_dims(ref_ref_model)
num_cell_dims(ref_ref_model)
get_panel_ids(ref_ref_model)
writevtk(get_ambient_grid(get_grid(ref_ref_model)),dir*"/ref_ref_grid",append=false)
