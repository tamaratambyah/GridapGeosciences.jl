"""
Set up the ThermogeostrophicBalanceTest from Section 5.1 of DOI. 10.1137/24M1638938
for TSW on the sphere

Note, I think equation (5.3) in the paper should have H^2
"""

using DrWatson
using JLD2

include("../convergence_tools.jl")
include("Williamson_functions.jl")


a_e = 6.37e6 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 2.94e4/g #m
u_0 = 2*π*a_e/T #m/s
c = 0.05
TF = 5*(3600*24)

L = a_e
_τ = 1/ω #sqrt(a_e/g)

_a = a_e/L
_g = g*_τ^2/L
_ω = ω*_τ
_H_0 = H_0/L
_T = T/_τ
_u0 = u_0/L*_τ #2*π*_a/_T
_tF = TF/_τ


function b₀(ζ)
  function _b₀(xyz)
    h = h₀(ζ)(xyz)
    _g*( 1 + c *_H_0^2 /h^2)
  end
end

function B₀(ζ)
  function _B₀(xyz)
    h = h₀(ζ)(xyz)
    b = b₀(ζ)(xyz)
    b*h
  end
end
