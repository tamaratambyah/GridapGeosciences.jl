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

"""
Note: the indexing of TensorValues is backwards to matrices. Index the "column"
first, and then the "row". A 2x3 TensorValue has the from
v = [1 2
     3 4
     5 6]
where
v[1,1] = 1;   v[2,1] = 2
v[1,2] = 3;   v[2,2] = 4
v[1,3] = 5;   v[2,3] = 6
"""
"""
_unit_outward_normal -- takes the cross product of the basis vectors of the
tangent space, which is perpendicular to the tangent space. i.e. normal to the
surface
"""
function _unit_outward_normal(v::Gridap.TensorValues.MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  n/norm(n)
end

### Analytical map
cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
CSgrid_analytical_H = CubedSphereGrid(cube_grid,map_cube_to_sphere)
CSmodel = CubedSphereDiscreteModel(CSgrid_analytical_H,topo,face_labels)

cell_reffe = get_reffes(get_grid(CSmodel))
cell_dofs = lazy_map(get_dof_basis,cell_reffe)
cell_x = lazy_map(get_nodes,cell_dofs)
cell_map = get_cell_map(CSmodel)
cell_Jt = lazy_map(Broadcasting(∇),cell_map)

_v = map(evaluate,cell_Jt[2],cell_x)
v = _v[1][1]

n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]

# Jt = lazy_map(∇,_cell_map)
n = lazy_map( Operation(_unit_outward_normal), cell_Jt)

map(n) do n
  println( map(evaluate,n,cell_x) )
end
nx = map(evaluate,n[3],cell_x)

normals = Gridap.CellData.GenericCellField(n,Triangulation(CSmodel),ReferenceDomain())

writevtk(Triangulation(CSmodel),dir*"/normals",cellfields=["n"=>normals],append=false)


"""
The columns of the Jacobian form a basis of the tangent space.
The tangent space is a vector space of dimension 2
"""
function get_tangent_space_basis(model::CubedSphereDiscreteModel)
  cell_reffe = get_reffes(get_grid(model))
  cell_dofs = lazy_map(get_dof_basis,cell_reffe)
  cell_x = lazy_map(get_nodes,cell_dofs)
  cell_map = get_cell_map(model)
  cell_Jt = lazy_map(Broadcasting(∇),cell_map)
  return cell_Jt, cell_x
end

function get_outward_normal_vector(trian::Triangulation)

  n = lazy_map( Operation(_unit_outward_normal), cell_Jt)

end
