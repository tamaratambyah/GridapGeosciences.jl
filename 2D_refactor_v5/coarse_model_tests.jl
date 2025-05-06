using Gridap
include("initalise.jl")

manifold_model = ManifoldDiscreteModel(cube_model_3D,cube)

num_point_dims(manifold_model)
num_cell_dims(manifold_model)

manifold_grid = get_grid(manifold_model)
topo = get_grid_topology(manifold_model)

get_cell_permutations(topo,0)
get_cell_node_ids(manifold_grid)

test_cell_maps(get_cell_map(manifold_grid),get_cell_ref_coordinates(manifold_grid),get_cell_coordinates(manifold_grid))

writevtk(manifold_grid,dir*"/grid",append=false)

panel_ids = get_panel_ids(manifold_model)
ambient_model = get_ambient_model(manifold_model)
test_cell_maps(ambient_model)

Ω_ambient = Triangulation(ambient_model)

function alpha(p)
  function _alpha(x)
    R = rp1_3D[p]
    y = R⋅x
    z = A_bump ⋅y
    z[1]
  end
end

function beta(p)
  function _beta(x)
    R = rp1_3D[p]
    y = R⋅x
    z = A_bump ⋅y
    z[2]
  end
end

pts = get_cell_points(Ω_ambient)

cell_field_a = map(p->GenericField(alpha(p)),panel_ids)
cf_alpha = GenericCellField(cell_field_a,Ω_ambient,PhysicalDomain())
cvals_alpha = cf_alpha(pts)


cell_field_b = map(p->GenericField(beta(p)),panel_ids)
cf_beta = GenericCellField(cell_field_b,Ω_ambient,PhysicalDomain())
cvals_beta = cf_beta(pts)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["a"=>cf_alpha,"b"=>cf_beta],append=false)
