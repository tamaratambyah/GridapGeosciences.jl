using Gridap
using FillArrays

function setup_model()
    ptr     = Int32[ 1, 5, 9 ]
    data    = Int32[ 1, 2, 3, 4, 2, 5, 4, 6 ]

    node_coordinates = Vector{Point{2,Float64}}(undef,6)


    cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
    cell_vertex_coordinates = Vector{Point{2,Float64}}(undef,8)

    cell_vertex_coordinates[1]=Point{2,Float64}(0.0,0.0)
    cell_vertex_coordinates[2]=Point{2,Float64}(1.0,0.0)
    cell_vertex_coordinates[3]=Point{2,Float64}(0.0,1.0)
    cell_vertex_coordinates[4]=Point{2,Float64}(1.0,1.0)

    cell_vertex_coordinates[5]=Point{2,Float64}(1.0,1.0)
    cell_vertex_coordinates[6]=Point{2,Float64}(1.0,0.0)
    cell_vertex_coordinates[7]=Point{2,Float64}(0.0,1.0)
    cell_vertex_coordinates[8]=Point{2,Float64}(0.0,0.0)

    cell_vertex_coordinates = Gridap.Arrays.Table(cell_vertex_coordinates, ptr)


    polytope=QUAD
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,2))
    cell_reffes=[scalar_reffe]

    cell_shape_funs = FillArrays.Fill( Gridap.ReferenceFEs.get_shapefuns(scalar_reffe), 2 )

    cell_map = lazy_map(Gridap.ReferenceFEs.linear_combination,cell_vertex_coordinates,cell_shape_funs)

    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.NonOriented(),
                                            nothing,
                                            cell_map)
    Gridap.Geometry.UnstructuredDiscreteModel(grid)
end

model = setup_model()
Ω = Triangulation(model)
dΩ = Measure(Ω,2)


coords = [Point(2.0,1.0),Point(2.0,0.0),Point(1.0,1.0),Point(1.0,0.0)]
polytope=QUAD
scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
shfuns = Gridap.ReferenceFEs.get_shapefuns(scalar_reffe)
mapx = Gridap.ReferenceFEs.linear_combination(coords,shfuns)

function f(p,x)
  if p==1
      return VectorValue(x[1],x[2])
  else
      x=mapx(x)
      return VectorValue(-x[2],-x[1])
  end
end

function g(p)
  x->f(p,x)
end

cf = Gridap.CellData.GenericCellField(map(p->Gridap.Fields.GenericField(g(p)),[1,2]),Ω, PhysicalDomain())

p_fe = 1
R = TestFESpace(model, ReferenceFE(lagrangian,VectorValue{2,Float64},p_fe);conformity=:L2)
get_cell_dof_ids(R)

#### Build FE space with constraints ###

# Dofs sitting on the first cell,
# lower right corner (2,6) and
# right upper corner (4,8) are the "slave" DOFs
sDOF_to_dof = [2,6,4,8]
# Dofs sitting on the second cell,
# lower left corner (9,13) and
# upper left corner (11,15)
# are the "master" DoFs. We swap and
# flip the sign of the master DOFs to get the slave DOFs
# (this is the change of basis among coordinate systems of the two cells)
sDOF_to_dofs = Gridap.Arrays.Table([ [13], [9], [15], [11] ])
sDOF_to_coeffs = Gridap.Arrays.Table([ [-1.0], [-1.0], [-1.0], [-1.0] ])
Rconstrained = Gridap.FESpaces.FESpaceWithLinearConstraints(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, R)

fh = interpolate(cf, Rconstrained)

eh=fh-cf
dc = ∫(eh⋅eh)*dΩ
