function xyz2θϕr(x)
  r = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
  θ = atan(x[2], x[1])
  ϕ = asin(x[3]/r)
  VectorValue(θ,ϕ,r)
end

function θϕ2xyz(θϕ)
  θ,ϕ = θϕ
  x = cos(θ)*cos(ϕ)
  y = sin(θ)*cos(ϕ)
  z = sin(ϕ)
  VectorValue(x,y,z)
end

function spherical_to_cartesian_matrix(θϕr)
  θ,ϕ,r = θϕr
  TensorValue(-sin(θ)       , cos(θ)       ,      0,
              -sin(ϕ)*cos(θ),-sin(ϕ)*sin(θ), cos(ϕ),
               cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), sin(ϕ))
end

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

# Topography
function topography(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  #θc    = -π/2.0
  θc    =  0.0
  ϕc    =  π/6.0
  rad   = π/9.0
  rsq   = (ϕ - ϕc)*(ϕ - ϕc) + (θ - θc)*(θ - θc)
  r     = sqrt(rsq)
  b     = 0.0
  if(r < rad)
    b = _b0*(1.0 - r/rad)
  end
  b
end
function _topography(xyz)
  0.0
end


a_e = 6.37e6 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 2.94e4/g #m
u_0 = 20 #2*π*a_e/T #m/s
b_0 = 2000.0 #m
TF = 1*(3600*24)

L = a_e
_τ = 1/ω#sqrt(a_e/g)

_a = a_e/L
_g = g*_τ^2/L
_ω = ω*_τ
_H_0 = H_0/L
_T = T/_τ
_u0 = u_0/L*_τ #2*π*_a/_T
_b0 = b_0/L
_tF = TF/_τ
