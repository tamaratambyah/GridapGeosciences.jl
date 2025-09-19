### Williamson2 functions
vWilliamson(ζ,u0,ω,grav,H0) = θϕ -> u0*VectorValue( cos(θϕ[2])*cos(ζ) + cos(θϕ[1])*sin(θϕ[2])*sin(ζ),
                                            -sin(θϕ[1])*sin(ζ) )

hWilliamson(ζ,u0,ω,grav,H0) = θϕ -> H0 - 1/grav*(ω*u0 + 0.5*u0^2)*( -cos(θϕ[1])*cos(θϕ[2])*sin(ζ) +  sin(θϕ[2])*cos(ζ) )^2

fWilliamson(ζ,u0,ω,grav,H0) = θϕ -> 2*ω*( - cos(θϕ[1])*cos(θϕ[2])*sin(ζ) + sin(θϕ[2])*cos(ζ) )


a_e = 6.37e6
g = 9.8
ω = 7.29e-5
T = 12*24*3600
H_0 = 2.94e4/g
u_0 = 2*π*a_e/T

L = a_e
τ = sqrt(a_e/g)

_a = 1.0
_g = 1.0
_ω = ω*τ
_H_0 = H_0/L
_T = T/τ
_u0 = 2*π*_a/_T

_vWilliamson(ζ) = θϕ -> _u0*VectorValue( cos(θϕ[2])*cos(ζ) + cos(θϕ[1])*sin(θϕ[2])*sin(ζ),
                                            -sin(θϕ[1])*sin(ζ) )

_hWilliamson(ζ) = θϕ -> _H_0 - (2*π/_T)*(_ω + π/_T)*( -cos(θϕ[1])*cos(θϕ[2])*sin(ζ) +  sin(θϕ[2])*cos(ζ) )^2

_fWilliamson(ζ) = θϕ -> 2*_ω*( - cos(θϕ[1])*cos(θϕ[2])*sin(ζ) + sin(θϕ[2])*cos(ζ) )

_ηWilliamson(ζ) = θϕ -> ( 4*π/_T + 2*_ω  )*( -cos(θϕ[1])*cos(θϕ[2])*sin(ζ) + sin(θϕ[2])*cos(ζ)  )

# _qWilliamson(ζ) = θϕ -> _ηWilliamson(ζ)(θϕ)/_hWilliamson(ζ)(θϕ)
_qWilliamson(ζ) = θϕ -> ( 4*π/_T + 2*_ω  )*( -cos(θϕ[1])*cos(θϕ[2])*sin(ζ) + sin(θϕ[2])*cos(ζ)  )/(_H_0 - (2*π/_T)*(_ω + π/_T)*( -cos(θϕ[1])*cos(θϕ[2])*sin(ζ) +  sin(θϕ[2])*cos(ζ) )^2)
