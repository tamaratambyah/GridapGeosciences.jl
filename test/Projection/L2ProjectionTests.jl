using Gridap
using Gridap.Arrays
using Gridap.Helpers
using GridapGeosciences
using Test

include("L2_projection_Hcurl.jl")
include("L2_projection_Hdiv.jl")
include("L2_projection_Lagrangian_scalar.jl")
include("L2_projection_Lagrangian_vector.jl")

# Vector field in the tangent space of the sphere
function uX(p)
  function _u(α)
    x = ForwardMap(p)(α)
    VectorValue(-x[2],x[1],0.0)
  end
end

# Scalar function
function fS(p)
  function f(αβ)
    xyz = ForwardMap(p)(αβ)
    xyz[1]*xyz[2]*xyz[3]
  end
end

function L2_projection(models::AbstractArray; _i_am_main=true)
  ls = LUSolver()
  dir = @__DIR__
  ps = [1,2]
  Dc = num_cell_dims(testitem(models))

  # L2: scalar
  p_convergence_auto_test(ps,models,L2_projection_Lagrangian_scalar,dir,fS,:L2,ls,false; _i_am_main=_i_am_main)

  # H1: scalar
  p_convergence_auto_test(ps,models,L2_projection_Lagrangian_scalar,dir,fS,:H1,ls,false; _i_am_main=_i_am_main)

  # L2: vector
  p_convergence_auto_test(ps,models,L2_projection_Lagrangian_vector,dir,uX,:L2,ls,false; _i_am_main=_i_am_main)

  # H1: vector
  p_convergence_auto_test(ps,models,L2_projection_Lagrangian_vector,dir,uX,:H1,ls,false; _i_am_main=_i_am_main)

  # H div: vector
  p_convergence_auto_test(ps,models,L2_projection_Hdiv,dir,uX,ls,false; _i_am_main=_i_am_main)

  # H curl: vector (3D)
  if Dc == 3
    ps = [0,1]
    p_convergence_auto_test(ps,models,L2_projection_Hcurl,dir,uX,ls,false; _i_am_main=_i_am_main)
  end
  @test true
end
