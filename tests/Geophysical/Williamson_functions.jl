# Coriolis
function f₀(ζ)
  function _f₀(xyz)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    2.0*_ω*( -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ) )
  end
end

# Initial velocity
function u₀(ζ)
  function _u₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  u     = _u0*(cos(ϕ)*cos(ζ) + cos(θ)*sin(ϕ)*sin(ζ))
  v     = - _u0*sin(θ)*sin(ζ)
  spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,v,0)
  end
end


# Initial fluid depth
function h₀(ζ)
  function _h₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  h  = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
  _H_0 - (_ω*_u0 + 0.5*_u0*_u0)*h*h/_g
  end
end

function η₀(ζ)
  function _η₀(xyz)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    η = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
    ( 4*π/_T + 2*_ω  )* η
  end
end

# q₀(xyz) = η₀(xyz)/ h₀(xyz)
function q₀(xyz) ### do not use
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  η = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
  h  = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)

  top = ( 4*π/_T + 2*_ω  )* η
  bottom = _H_0 - (_ω*_u0 + 0.5*_u0*_u0)*h*h/_g

  top/bottom
end
