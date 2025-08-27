### Williamson2 functions
vWilliamson(ζ,u0,ω) = θϕ -> u0*VectorValue( cos(θϕ[2])*cos(ζ) + cos(θϕ[1])*sin(θϕ[2])*sin(ζ),
                                            -sin(θϕ[1])*sin(ζ) )

hWilliamson(ζ,u0,ω) = θϕ -> 1 - (ω*u0 + 0.5*u0^2)*( -cos(θϕ[1])*cos(θϕ[2])*sin(ζ) +  sin(θϕ[2])*cos(ζ) )^2

fWilliamson(ζ,u0,ω) = θϕ -> 2*ω*( - cos(θϕ[1])*cos(θϕ[2])*sin(ζ) + sin(θϕ[2])*cos(ζ) )
