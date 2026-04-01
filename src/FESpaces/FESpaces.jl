module FESpaces

# BEGIN of copy/paste from Gridap 0.20.x
# This part has been copy/pasted from Gridap 0.20.x. 
# Once we move to Gridap 0.20 from 0.19, we can remove this part of the 
# module and use it from Gridap.
import Gridap.Fields: LinearCombinationMap
import Gridap.Fields: linear_combination
import Gridap.Arrays: return_cache, evaluate, evaluate!
import Gridap.FESpaces: Dof
import Gridap.Fields: Point
import Gridap.ReferenceFEs: PointValue
import Gridap.Helpers: @check
include("LinearCombinationDofVector.jl")
# END of copy/paste from Gridap 0.20.x

import Gridap.ReferenceFEs: GenericLagrangianRefFE, GradConformity, get_prebasis, 
                            get_shapefuns, get_node_coordinates, get_polytope, num_dofs,
                            get_face_own_dofs, get_face_own_nodes, num_faces, get_offset, 
                            get_dof_basis

import Gridap.Geometry: get_grid_topology, num_faces, num_cell_dims, 
                        get_faces, num_cells, get_cell_map, 
                        BodyFittedTriangulation
import Gridap.Fields: Map, VectorValue, TensorValue, Point
import Gridap.Arrays: array_cache, CachedMatrix, return_type, getindex!, testitem, lazy_map, IdentityVector, setsize!
import GridapGeosciences.Geometry: ParametricDiscreteModel, get_panel_ids
import GridapGeosciences.Helpers: forward_jacobian, forward_pinv_jacobian
import Gridap.FESpaces: get_cell_shapefuns, get_cell_dof_basis, _use_clagrangian, H1Conformity
import Gridap.Helpers: @notimplemented
import LinearAlgebra: I, ⋅
include("GradConformingFESpaces.jl")
export _generate_change_of_basis_matrices

end