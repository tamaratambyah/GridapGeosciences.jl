# Solves the Galewsky test case for the shallow water equations on a sphere
# of physical radius 6371220m. Involves a shear flow instability of a zonal
# jet triggered by an initial gravity wave.
# reference:
#   Galewsky, Scott and Polvani (2004) Tellus, 56A 429-440

using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapGeosciences
using GridapGeosciences.Distributed


a_e = 6.37e6 # m
rₑ = a_e #m
g = 9.8 # m/2
Ωₑ  = 7.292e-5 #1/s
H₀ = 10000.0 #m
h0 = 120.0 #m
umax  = 80.0 #m/s
TF = 20*(3600*24) #s

L = a_e
_τ = 1/Ωₑ

_a = a_e/L
_rₑ = rₑ/L
_g = g*_τ^2/L
_Ωₑ = Ωₑ*_τ
_H_0 = H₀/L
_h0 = h0/L
_umax  = umax /L*_τ
_tF = TF/_τ




function uθ(θϕr)
  θ,ϕ,r = θϕr
  ϵ     = 1.0e-8
  ϕ₁    = π/7
  ϕ₂    = π/2 - ϕ₁
  en    = exp(-4.0/((ϕ₂ - ϕ₁)*(ϕ₂ - ϕ₁)))
  u     = 0.0
  if ϕ > ϕ₁ + ϵ && ϕ < ϕ₂ - ϵ
    u = (_umax/en)*exp(1.0/((ϕ - ϕ₁)*(ϕ - ϕ₂)))
  end
  u
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
  h     = _H_0
  hh    = _h0
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
    _f   = 2.0*_Ωₑ*sin(ϕₚ)
    h    = h - _rₑ*u*(_f + tan(ϕₚ)*u/_rₑ)*dϕ/_g
  end
  h = h + hh*cos(ϕ)*exp(-1.0*(θ/α)*(θ/α))*exp(-1.0*((ϕ₂ - ϕ)/β)*((ϕ₂ - ϕ)/β))
  h
end


function f₀(x)         # Coriolis term
  2.0*_Ωₑ*x[3]/_rₑ
end


function b₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  _g*(1 - 0.1*cos(ϕ)*exp( -(3*θ)^2 - (15 * (π/4-ϕ) )^2  )   )
end

function B₀(xyz)
  h = h₀(xyz)
  b = b₀(xyz)
  h*b
end





# dir = datadir("Galewsky")
# !isdir(dir) && mkdir(dir)

# MPI.Init()
# ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

# n_ref_lvls = 5

# o3model = GridapGeosciences.Distributed.ParametricOctreeDistributedDiscreteModel(ranks;
#   num_initial_uniform_refinements=n_ref_lvls)
# panel_model = o3model.parametric_dmodel


# panel_ids = get_panel_ids(panel_model)
# Ω_panel = Triangulation(panel_model)

# vX = panel_to_cartesian(tangent_vec(u₀))
# f = panel_to_cartesian(f₀)
# h = panel_to_cartesian(h₀)
# B = panel_to_cartesian(B₀)

# u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
# f_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
# h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
# B_cf = panelwise_cellfield(B,Ω_panel,panel_ids)
# covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)


# writevtk(Ω_panel,dir*"/IC",
#   cellfields=["u"=>covarient_basis_cf⋅u_contra_cf, "f"=>f_cf, "h"=>h_cf, "B"=>B_cf],
# append=false, geo_map=geo_map_func(Ω_panel))
