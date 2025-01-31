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
include("cube_surface_1_cell_per_panel.jl")
include("structs.jl")
include("refinement.jl")
include("maps.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()



### Polynomial map
CSgrid_poly = CubedSphereGrid(cube_grid,map_cube_to_latlon,1)
CSmodel_poly = CubedSphereDiscreteModel(CSgrid_poly,topo,face_labels)
writevtk(CSmodel_poly,dir*"/CSmodel_poly",append=false)

# test node 4 = (1,1,1) is mapped properly in cells 1,2,3 for mapping of order 1
cmaps_poly = collect( get_cell_map(CSmodel_poly) )
@test evaluate(cmaps_poly[1],Point(1,1)) == mapped_node4
@test evaluate(cmaps_poly[2],Point(1,0)) == mapped_node4
@test evaluate(cmaps_poly[3],Point(0,0)) == mapped_node4

# apply order 2 mapping
CSgrid_poly2 = CubedSphereGrid(cube_grid,map_cube_to_latlon,2,transfer=true)
CSmodel_poly2 = CubedSphereDiscreteModel(CSgrid_poly2,topo,face_labels)

# apply refinement
CSmodel_poly2_refined1 = Gridap.Adaptivity.refine(CSmodel_poly2)
writevtk(CSmodel_poly2_refined1,dir*"/CSmodel_poly_refined1",append=false)

CSmodel_poly2_refined2 = Gridap.Adaptivity.refine(CSmodel_poly2_refined1)
writevtk(CSmodel_poly2_refined2,dir*"/CSmodel_poly_refined2",append=false)


model = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
order = 1
 # Create a FESpace for the geometrical description
  # The FEFunction that describes the coordinate field
  # or displacement can be a interpolation or a solution
  # of a mesh displacement problem
  T = eltype(get_node_coordinates(model))
  Ts = eltype(T)

  os = orders
  _ps = get_polytopes(model)
  ct = get_cell_type(model)
  ps = lazy_map(Reindex(_ps), ct)


  Vₕ = FESpace(model,ReferenceFE(lagrangian,T,order);conformity=:H1)


  Vₕ_scal = FESpace(model,ReferenceFE(lagrangian,Ts,order);conformity=:H1)

  grid = get_grid(model)
  geo_map = get_cell_map(grid)

  cell_ctype = get_cell_type(grid)
  c_reffes = get_reffes(grid)

  # Create a fe_map using the cell_map that can be evaluated at the
  # vertices of the fe space (nodal type)
  # This returns a FEFunction initialised with the coordinates
  # But this is to build the FEFunction that will be inserted, it is
  # an advanced constructor, not needed at this stage

  c_dofs = get_fe_dof_basis(Vₕ)
  dof_basis = get_data(c_dofs)
  c_nodes = lazy_map(get_nodes,get_data(c_dofs))

  xh = zero(Vₕ)
  c_dofv = lazy_map(evaluate,dof_basis,geo_map)

  Uₕ = TrialFESpace(Vₕ)

  fv = get_free_dof_values(xh)
  dv = get_dirichlet_dof_values(get_fe_space(xh))
  gather_free_and_dirichlet_values!(fv,dv,Uₕ,c_dofv)

  c_xh = lazy_map(evaluate,get_data(xh),c_nodes)
  c_scal_ids = get_cell_dof_ids(Vₕ_scal)


  nodes_coords = Vector{eltype(eltype(c_xh))}(undef,num_free_dofs(Vₕ_scal))
  Geometry._cell_vector_to_dof_vector!(nodes_coords,c_scal_ids,c_xh)

  grid_map = GridWithFEMap(grid, Vₕ, Vₕ_scal, xh, nodes_coords, reffes)

_Vₕ = FESpace(model,ReferenceFE(lagrangian,VectorValue{2,Float64},order);conformity=:H1)

FE_map = interpolate(map_cube_to_latlon,  _Vₕ)


  add_mesh_displacement!(grid_map,FE_map)
