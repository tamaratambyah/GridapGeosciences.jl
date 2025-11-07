# Initial velocity
function u₀(ζ)
  function _u₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  u     = _u0*(cos(ϕ)*cos(ζ) + cos(θ)*sin(ϕ)*sin(ζ))
  v     = - _u0*sin(θ)*sin(ζ)
  ## in P6est, the radial component is first
  spherical_to_cartesian_matrix_3D(θϕr)⋅VectorValue(0.0,u,v)
  end
end

"""
X = r cosθ cosϕ
Y = r sinθ sinϕ
Z = r sinϕ
J = [dXdr dXdθ dXdϕ
     dYdr dYdθ dYdϕ
     dZdr dZdθ dZdϕ ]
As a TensorValue:
TensorValue = (dXdr,dYdr,dZdr,  dXdθ,dYdθ,dZdθ,  dXdϕ,dYdϕ,dZdϕ)
"""
function spherical_to_cartesian_matrix_3D(θϕr)
  θ,ϕ,r = θϕr
  TensorValue(cos(θ)*cos(ϕ),     sin(θ)*cos(ϕ), sin(ϕ),
             -r*sin(θ)*cos(ϕ), r*cos(θ)*cos(ϕ), 0,
             -r*cos(θ)*sin(ϕ), -r*sin(θ)*sin(ϕ), r*cos(ϕ) )
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

# Coriolis
function f₀(ζ)
  function _f₀(xyz)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    2.0*_ω*( -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ) )
  end
end

a_e = 6.37e6 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 2.94e4/g #m
u_0 = 2*π*a_e/T #m/s
b_0 = 2000.0 #m
TF = 1*(3600*24)

L = a_e
_τ = sqrt(a_e/g)

_a = a_e/L
_g = g*_τ^2/L
_ω = ω*_τ
_H_0 = H_0/L
_T = T/_τ
_u0 = u_0/L*_τ #2*π*_a/_T
_b0 = b_0/L
_tF = TF/_τ
