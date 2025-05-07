using Gridap
include("../src/initialise.jl")

sphere_manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
num_point_dims(sphere_manifold_model)
num_cell_dims(sphere_manifold_model)

manifold_grid = get_grid(sphere_manifold_model)
topo = get_grid_topology(sphere_manifold_model)

get_cell_permutations(topo,0) # check no permuations
get_cell_node_ids(manifold_grid)

test_cell_maps(get_cell_map(manifold_grid),get_cell_ref_coordinates(manifold_grid),get_cell_coordinates(manifold_grid))

writevtk(manifold_grid,dir*"/grid",append=false)

################################################################################
#### test the ambient grid is the same as a cube with faces [-1,1]
################################################################################

ambient_grid = get_ambient_grid(manifold_grid)
ambient_cell_coords = get_cell_coordinates(ambient_grid)

round_ambient_cell_coords = lazy_map(RoundVectorValues(),ambient_cell_coords)

_cube_grid_3D = _coarse_cube_model_3D(1.0)
cube_cell_coords = get_cell_coordinates(_cube_grid_3D)

cube_cell_coords == round_ambient_cell_coords


################################################################################
#### plot alpha,beta in the ambient space to check the orientation
################################################################################
sphere_manifold_model = Adaptivity.refine(sphere_manifold_model)

panel_ids = get_panel_ids(sphere_manifold_model)
ambient_model = get_ambient_model(sphere_manifold_model)
test_cell_maps(ambient_model)

Ω_ambient = Triangulation(ambient_model)
pts = get_cell_points(Ω_ambient)

function alpha(p)
  function _alpha(x)
    R = rp1_3D[p]
    X = R⋅x # Xs,Ys,Zs on panel 1
    latlon = evaluate(Sigma(),X) # latlon on panel 1
    z = evaluate(InvGnomonicMap(),latlon)
    z[1]
  end
end

function beta(p)
  function _beta(x)
    R = rp1_3D[p]
    X = R⋅x # Xs,Ys,Zs on panel 1
    latlon = evaluate(Sigma(),X) # latlon on panel 1
    z = evaluate(InvGnomonicMap(),latlon)
    z[2]
  end
end


function lon(x)
  X,Y,Z = x
  # acos(X/sqrt(X^2+Y^2))
  # asin(Y/sqrt(X^2+Y^2))
  # atan(Y,X)
  rem2pi(Float64(my_atan(X,Y)), RoundDown)
end

function lat(x)
  latlon = evaluate(Sigma(),x) # latlon on panel p
  latlon[2]
end

X(x) = x[1]
Y(x) = x[2]
Z(x) = x[3]

cell_field_a = map(p->GenericField(alpha(p)),panel_ids)
cf_alpha = GenericCellField(cell_field_a,Ω_ambient,PhysicalDomain())
cvals_alpha = cf_alpha(pts)

cell_field_b = map(p->GenericField(beta(p)),panel_ids)
cf_beta = GenericCellField(cell_field_b,Ω_ambient,PhysicalDomain())
cvals_beta = cf_beta(pts)

cf_lon = CellField(lon,Ω_ambient)
cvals_lon = cf_lon(pts)

cf_lat = CellField(lat,Ω_ambient)
cvals_lat = cf_lat(pts)

cf_X = CellField(X,Ω_ambient)
cf_Y = CellField(Y,Ω_ambient)
cf_Z = CellField(Z,Ω_ambient)

writevtk(Ω_ambient,dir*"/ambient",
        cellfields=["a"=>cf_alpha,"b"=>cf_beta, "lat"=>cf_lat, "lon"=>cf_lon,
        "X"=>cf_X,"Y"=>cf_Y,"Z"=>cf_Z],
        append=false)

################################################################################
#### create a Raviart Thomas space
################################################################################
V = TestFESpace(sphere_manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
f(x) = VectorValue(x[1],x[2])
uh = interpolate(f ,V)

Ω_parametric = Triangulation(sphere_manifold_model)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)
