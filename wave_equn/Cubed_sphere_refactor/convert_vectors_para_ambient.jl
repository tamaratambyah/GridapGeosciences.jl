using DrWatson
using Gridap
using GridapGeosciences
using Test



θϕ = VectorValue(π,π/4) # xyz2θϕ(x) # point in θϕ space
Gθϕ = GridapGeosciences.G_unit_sphere(θϕ)
Jθϕ = GridapGeosciences.J_unit_sphere(θϕ)
JT = transpose(Jθϕ)
x = θϕ2xyz(θϕ) #Point(1,1,1) # point on sphere

Jθϕ⋅θϕ
π/4*1/sqrt(2)
-π/sqrt(2)

JT⋅(Jθϕ⋅θϕ)
π/2
π/4


vθϕ(θϕ) = VectorValue(θϕ[1],θϕ[2])
vxyz(x) = VectorValue(x[1],x[2],x[3])

function parametric_2_ambient(vθϕ)
  function tmp(θϕ::Point{2})
    Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
    Jθϕ⋅(vθϕ)(θϕ)
  end
end

function ambient_2_parametric(vxyz)
  function tmp(x::Point{3})
    θϕ = xyz2θϕ(x)
    Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
    transpose(Jθϕ)⋅(vxyz)(x)
  end
end

function _parametric_2_ambient(vθϕ::VectorValue{2,Float64},θϕ::Point{2})
    Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
    Jθϕ⋅(vθϕ)
end

function _ambient_2_parametric(vxyz::VectorValue{3,Float64},x::Point{3})
    θϕ = xyz2θϕ(x)
    Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
    transpose(Jθϕ)⋅(vxyz)
end

evaluate(parametric_2_ambient(vθϕ),(θϕ))
evaluate(ambient_2_parametric(vxyz),(x) )
parametric_2_ambient(vθϕ)(θϕ)
ambient_2_parametric(vxyz)(x)

_parametric_2_ambient(vθϕ(θϕ),θϕ)
_ambient_2_parametric(vxyz(x),x)
