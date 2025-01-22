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

dir = datadir("CubedSphereRefactor/Normals")
!isdir(dir) && mkdir(dir)

cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()

### Analytical map
CSgrid_analytical_H = CubedSphereGrid(cube_grid,map_cube_to_cube)
CSmodel = CubedSphereDiscreteModel(CSgrid_analytical_H,topo,face_labels)


# model = UnstructuredDiscreteModel( CartesianDiscreteModel((-1,1,-1,1,-1,1),(1,1,1)) )
model = CSmodel

D = 2#num_cell_dims(model)-1# = 2
nfaces = num_faces(model,D)

topo = get_grid_topology(model)

face_to_mask = fill(true,nfaces) ##get_isboundary_face(topo,D-1)
face_to_bgface = findall(face_to_mask)
bgface_to_lcell = Fill(1,num_facets(model))

bgface_grid = Grid(ReferenceFE{D},model)

face_grid = view(bgface_grid,face_to_bgface)
cell_grid = get_grid(model)

glue = Gridap.Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)

trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)
Γ =  BoundaryTriangulation(trian,glue)



# function Gridap.Geometry.FaceToCellGlue(
#   topo::GridTopology,
#   cell_grid::CubedSphereGrid,
#   face_grid::Grid,
#   face_to_bgface::AbstractVector,
#   bgface_to_lcell::AbstractVector)
  println("my func")
  D = num_cell_dims(cell_grid)
  cD = num_cell_dims(face_grid)
  bgface_to_cells = get_faces(topo,cD,D)
  cell_to_bgfaces = get_faces(topo,D,cD)
  cell_to_lface_to_pindex = Table(get_cell_permutations(topo,cD))

  bgface_to_cell = lazy_map(getindex,bgface_to_cells, bgface_to_lcell)
  bgface_to_lface = find_local_index(bgface_to_cell, cell_to_bgfaces)

  face_to_cell = collect(Int32,lazy_map(Reindex(bgface_to_cell), face_to_bgface))
  face_to_lface = collect(Int8,lazy_map(Reindex(bgface_to_lface), face_to_bgface))
  face_to_lcell = collect(Int8,lazy_map(Reindex(bgface_to_lcell), face_to_bgface))

  f = (p)->fill(Int8(UNSET),num_faces(p,cD))
  ctype_to_lface_to_ftype = map( f, get_reffes(cell_grid) )
  face_to_ftype = get_cell_type(face_grid)
  cell_to_ctype = get_cell_type(cell_grid)

  _fill_ctype_to_lface_to_ftype!(
    ctype_to_lface_to_ftype,
    face_to_cell,
    face_to_lface,
    face_to_ftype,
    cell_to_ctype)

  FaceToCellGlue(
    face_to_bgface,
    bgface_to_lcell,
    face_to_cell,
    face_to_lface,
    face_to_lcell,
    face_to_ftype,
    cell_to_ctype,
    cell_to_lface_to_pindex,
    ctype_to_lface_to_ftype)
end

#### facet normal
trian = Γ
boundary_trian_glue = glue

cell_grid = get_grid(get_background_model(trian.trian))
_cell_grid = get_grid(CSmodel)


r =  get_reffes(cell_grid)[1]
p = Gridap.ReferenceFEs.get_polytope(r)
    lface_to_n = get_facet_normal(p)
    lface_to_pindex_to_perm = get_face_vertex_permutations(p,num_cell_dims(p)-1)
    nlfaces = length(lface_to_n)
    lface_pindex_to_n = [ fill(lface_to_n[lface],length(lface_to_pindex_to_perm[lface])) for lface in 1:nlfaces ]
    lface_pindex_to_n

  ## Reference normal
  function f(_r)
    r = LagrangianRefFE(Float64,HEX,1)
    p = get_polytope(r)
    lface_to_n = get_facet_normal(p)
    lface_to_pindex_to_perm = get_face_vertex_permutations(p,num_cell_dims(p)-1)
    nlfaces = length(lface_to_n)
    lface_pindex_to_n = [ fill(lface_to_n[lface],length(lface_to_pindex_to_perm[lface])) for lface in 1:nlfaces ]
    lface_pindex_to_n
  end
  ctype_lface_pindex_to_nref = map(f, get_reffes(cell_grid)) # same as cube
  face_to_nref = Gridap.Geometry.FaceCompressedVector(ctype_lface_pindex_to_nref,boundary_trian_glue)
  face_s_nref = lazy_map(Gridap.Fields.constant_field,face_to_nref)

  # Inverse of the Jacobian transpose
  cell_q_x = get_cell_map(cell_grid)
  cell_q_Jt = lazy_map(∇,cell_q_x)
  cell_q_invJt = lazy_map(Operation(Gridap.Fields.pinvJt),cell_q_Jt)
  face_q_invJt = lazy_map(Reindex(cell_q_invJt),boundary_trian_glue.face_to_cell)

  # Change of domain
  D = num_cell_dims(cell_grid)
  boundary_trian_glue = get_glue(trian,Val(D))
  face_s_q = boundary_trian_glue.tface_to_mface_map
  face_s_invJt = lazy_map(∘,face_q_invJt,face_s_q)
  face_s_n = lazy_map(Broadcasting(Operation(Gridap.Geometry.push_normal)),face_s_invJt,face_s_nref)
  Fields.MemoArray(face_s_n)
