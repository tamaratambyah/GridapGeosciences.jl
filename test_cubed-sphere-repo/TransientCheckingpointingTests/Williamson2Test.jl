using DrWatson

include(projectdir("test/Geophysical/Williamson_functions.jl"))

function topography(xyz)
  0.0
end


a_e = 6.37e6 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 2.94e4/g #m
u_0 = 2*π*a_e/T #m/s
b_0 = 0.0 #m
TF = 5*(3600*24)

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
