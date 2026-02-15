module ODEs
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces, Gridap.ODEs
using LinearAlgebra
using FillArrays
using PartitionedArrays

import PartitionedArrays: consistent!

import Gridap.Helpers: @abstractmethod
import Gridap.ODEs: LinearODE, QuasilinearODE, SemilinearODE
import Gridap.ODEs: ODEOperatorType
import Gridap.ODEs: ODEOperator

import Gridap.ODEs: Assembler
import Gridap.FESpaces: collect_cell_vector
import Gridap.FESpaces: collect_cell_matrix
import Gridap.FESpaces: assemble_vector_add!
import Gridap.FESpaces: allocate_matrix
import Gridap.FESpaces: assemble_matrix_add!
import Gridap.FESpaces: assemble_matrix!
import Gridap.FESpaces: assemble_vector!
import Gridap.ODEs: TransientFEOperator
import Gridap.Algebra: NonlinearSolver
import Gridap.FESpaces: FEOperatorFromWeakForm
import Gridap.FESpaces: FEOperator
import Gridap.ODEs: EXRungeKutta, DIMRungeKutta
import Gridap.ODEs: AbstractTableau
import Gridap.ODEs: AbstractQuasilinearODE, AbstractSemilinearODE
using LinearAlgebra
import Gridap.ODEs: get_weights, get_nodes,get_matrix
import Gridap.ODEs: RungeKutta
import Gridap.ODEs: ODESolver
import Gridap.ODEs: TableauType, ExplicitTableau, ImplicitTableau

import Gridap.ODEs: StageOperator
import Gridap.ODEs: LinearStageOperator, NonlinearStageOperator

include("DAE.jl")
include("ODEOperators.jl")
include("ODEOpsFromTFEOps.jl")
include("RungeKuttaEX.jl")
include("RungeKuttaDIM.jl")
include("StageOperators.jl")

export DAEFEOperator
export LinearStageOperator
export NonlinearStageOperator

end
