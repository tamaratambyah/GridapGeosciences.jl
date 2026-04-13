include("../Geophysical/Williamson_functions.jl")

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


a_e = 6.37e6 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 5960.0 #2.94e4/g #m
u_0 = 20 #2*π*a_e/T #m/s
b_0 = 2000.0 #m
TF = 20*(3600*24) # 15 days as per Rognes 2013 paper

L = a_e
_τ = 1/ω#

_a = a_e/L
_g = g*_τ^2/L
_ω = ω*_τ
_H_0 = H_0/L
_T = T/_τ
_u0 = u_0/L*_τ #2*π*_a/_T
_b0 = b_0/L
_tF = TF/_τ
