# Solves the Galewsky test case for the shallow water equations on a sphere
# of physical radius 6371220m. Involves a shear flow instability of a zonal
# jet triggered by an initial gravity wave.
# reference:
#   Galewsky, Scott and Polvani (2004) Tellus, 56A 429-440



a_e = 6.37e6 # m
rₑ = a_e #m
g = 9.8 # m/2
Ωₑ  = 7.292e-5 #1/s
H₀ = 10000.0 #m
h0 = 120.0 #m
umax  = 80.0 #m/s
# TF = 20*(3600*24) #s
TF = 0.1*(3600*24) #s

L = a_e
_τ = 1/Ωₑ
# _τ = sqrt(a_e/g)

Uscale = L/_τ # m/s
Hscale = L #m
gscale = L/_τ^2 # m/s^2
bscale = gscale # m/s^2
fscale = 1/_τ # 1/s

_g = g/gscale
_H_0 = H₀/Hscale
_tF = TF/_τ


function uθ(θϕr)
  θ,ϕ,r = θϕr
  ϵ     = 1.0e-8
  ϕ₁    = π/7
  ϕ₂    = π/2 - ϕ₁
  en    = exp(-4.0/((ϕ₂ - ϕ₁)*(ϕ₂ - ϕ₁)))
  u     = 0.0
  if ϕ > ϕ₁ + ϵ && ϕ < ϕ₂ - ϵ
    u = (umax/en)*exp(1.0/((ϕ - ϕ₁)*(ϕ - ϕ₂)))
  end
  u/Uscale
end

# Initial velocity
function u₀(xyz)
  θϕr = xyz2θϕr(xyz)
  u   = uθ(θϕr)
  _spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,0,0)
end

function _spherical_to_cartesian_matrix(θϕr)
  θ,ϕ,r = θϕr
  TensorValue(-sin(θ)       , cos(θ)       ,      0,
              -sin(ϕ)*cos(θ),-sin(ϕ)*sin(θ), cos(ϕ),
               cos(ϕ)*cos(θ), cos(ϕ)*sin(θ), sin(ϕ))
end

# Initial fluid depth
function h₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  x,y,z = xyz
  θ,ϕ,r = θϕr
  h     = H₀
  hh    = 120.0
  α     = 1.0/3.0
  β     = 1.0/15.0
  ϕ₂    = π/4
  ni    = 1000
  ϕₚ    = 0.0
  dϕ    = abs(ϕ/ni)
  sgn   = 1.0
  if ϕ < 0.0
    sgn = -1.0
  end
  for i in 1:ni
    ϕₚ   = ϕₚ + sgn*dϕ
    _θϕr = VectorValue(θ,ϕₚ,r)
    u    = uθ(_θϕr)
    _f   = 2.0*Ωₑ*sin(ϕₚ)
    h    = h - rₑ*u*(_f + tan(ϕₚ)*u/rₑ)*dϕ/g
  end
  h = h + hh*cos(ϕ)*exp(-1.0*(θ/α)*(θ/α))*exp(-1.0*((ϕ₂ - ϕ)/β)*((ϕ₂ - ϕ)/β))
  h/Hscale
end


function f₀(x)         # Coriolis term
  f = 2.0*Ωₑ*x[3]/rₑ
  f/fscale
end


function b₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  # b = g
  b = g*(1 - 0.1*cos(ϕ)*exp( -(3*θ)^2 - (15 * (π/4-ϕ) )^2  )   )
  b/gscale
end

function B₀(xyz)
  h = h₀(xyz)
  b = b₀(xyz)
  h*b
end
