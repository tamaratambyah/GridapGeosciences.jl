
# Initial fluid depth
function h_3D(p)
  function _h(γαβ)
    ζ = 0.0
    xyz = forward_map_3D(p)(γαβ)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    h  = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
    _H_0 - (_ω*_u0 + 0.5*_u0*_u0)*h*h/_g
  end
end

function f_3D(p)
  function _f(γαβ)
    ζ = 0.0
    xyz = forward_map_3D(p)(γαβ)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    2.0*_ω*( -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ) )
  end
end

function f_vec_3D(p)
  function _f(γαβ)
    f = f_3D(p)(γαβ)

    xyz = forward_map_3D(p)(γαβ)
    n = normal_vec(xyz)
    f*n
  end
end



function η_3D(p)
  function _η(γαβ)
    ζ = 0.0
    xyz = forward_map_3D(p)(γαβ)
    θϕr   = xyz2θϕr(xyz)
    θ,ϕ,r = θϕr
    η = -cos(θ)*cos(ϕ)*sin(ζ) + sin(ϕ)*cos(ζ)
    ( 4*π/_T + 2*_ω  )* η
  end
end


function η_vec_3D(p)
  function _η(γαβ)
    η = η_3D(p)(γαβ)

    xyz = forward_map_3D(p)(γαβ)
    n = normal_vec(xyz)
    η*n
  end
end


function u_vec_3D(p)
  function _u(γαβ)
    ζ = 0.0
    xyz = forward_map_3D(p)(γαβ)
    # θϕr   = xyz2θϕr(xyz)
    # θ,ϕ,r = θϕr
    # u     = _u0*(cos(ϕ)*cos(ζ) + cos(θ)*sin(ϕ)*sin(ζ))
    # v     = - _u0*sin(θ)*sin(ζ)
    # _spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,v,0)

    #### from Rognes2013 paper
    r = sqrt(xyz[1]^2 + xyz[2]^2 + xyz[3]^2)
    u = -_u0*xyz[2]/r
    v = _u0*xyz[1]/r
    w = 0.0
    VectorValue(u,v,w)
  end
end

function _spherical_to_cartesian_matrix(θϕr)
  θ,ϕ,r = θϕr
  TensorValue(-r*sin(θ)*cos(ϕ), r*cos(θ)*cos(ϕ),      0,
              -r*sin(ϕ)*cos(θ),-r*sin(ϕ)*sin(θ), r*cos(ϕ),
               cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), sin(ϕ))
end

# function _spherical_to_cartesian_matrix(θϕr)
#   θ,ϕ,r = θϕr
#   TensorValue(-sin(θ)       , cos(θ)       ,      0,
#               -sin(ϕ)*cos(θ),-sin(ϕ)*sin(θ), cos(ϕ),
#                cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), sin(ϕ))
# end

function topography(p)
  function _u(γαβ)
    xyz = forward_map_3D(p)(γαβ)
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
end

#### Parameters for a reduced earth
a_e = 6.37e6#/125 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 5960.0 #2.94e4/g #m
u_0 = 20 #2*π*a_e/T #m/s
b_0 = 2000.0 #m
TF = 20*(3600*24)

L = a_e
_τ = 1/ω #sqrt(a_e/g)

_a = a_e/L
_g = g*_τ^2/L
_ω = ω*_τ
_H_0 = H_0/L
_T = T/_τ
_u0 = u_0/L*_τ #2*π*_a/_T
_b0 = b_0/L
_tF = TF/_τ
