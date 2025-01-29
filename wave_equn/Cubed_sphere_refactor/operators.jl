using DrWatson
using Gridap
using GridapGeosciences
using Test


x = Point(1,1,1) # point on sphere
θϕ = xyz2θϕ(x) # point in θϕ space
Gθϕ = GridapGeosciences.G_unit_sphere(θϕ)
Jθϕ = GridapGeosciences.J_unit_sphere(θϕ)
JT = transpose(Jθϕ)



### scalar functions
uxyz(x::Point{3}) = exp(x[1]) + sin(x[2]^2) + x[3] # xyz
uxyz(x)
function uθϕ(θϕ::Point{2},u)
  xyz = θϕ2xyz(θϕ)
  u(xyz)
end
uθϕ(θϕ::Point{2}) = uθϕ(θϕ,uxyz)
u = uθϕ(θϕ)

### vector functions
wxyz(x::Point{3}) = VectorValue((x[1]),x[2]^3,x[3])
wxyz(x)
_wθϕ(θϕ::Point{2}) = uθϕ(θϕ,wxyz)
_wθϕ(θϕ)

wθϕ(θϕ) = w_xyz2θϕ(x,wxyz)(θϕ)
wθϕ(θϕ)

wx(x) = w_θϕ2xyz(θϕ,wθϕ)(x)
wx(x)

function w_xyz2θϕ(x::Point{3},wx)
  function tmp(θϕ)
    Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
    JT = transpose(Jθϕ)
    JT⋅wx(x)
  end
end

function w_θϕ2xyz(θϕ::Point{2},wθϕ)
  Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
  function tmp(x)
    Jθϕ ⋅wθϕ(x)
  end
end



### gradient
_uθϕ(θϕ) = θϕ[1] + θϕ[2]
∇u(θϕ) = gradient_unit_sphere(_uθϕ)(θϕ) # in ambient space
∇u(θϕ)

"""
* input parametric v
* output parametric grad v
"""
function parametric_gradient_unit_sphere(vθϕ)
  function tmp(θϕ)
     Gθϕ=GridapGeosciences.G_unit_sphere(θϕ)
     inv(Gθϕ)⋅(∇(vθϕ))(θϕ)
  end
end

"""
* input paramtric v
* output ambient grad v
"""
function ambient_gradient_unit_sphere(vθϕ)
  function tmp(θϕ)
    Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
    Jθϕ⋅parametric_gradient_unit_sphere(vθϕ)(θϕ)
  end
end



parametric_gradient_unit_sphere(_uθϕ)(θϕ) # in parametric space
ambient_gradient_unit_sphere(_uθϕ)(θϕ) # in ambient space

ambient(x) = ambient_vector(parametric_gradient_unit_sphere(_uθϕ) )(θϕ)
parametric_vector(ambient)(x)

#### fixing divergence##################################################

w(x) = VectorValue(1,1,1)
_wx = w(x)
f(θϕ) =   transpose(Jθϕ)⋅(w(θϕ))
f(θϕ)
# function f(θϕ)
#   sqrt(det(Gθϕ))*inv(Gθϕ)⋅transpose(Jθϕ)⋅(v(θϕ))
# end
1.0/sqrt(det(Gθϕ))*(∇⋅(f))(θϕ)

#########################################################3

"""
* input paramtric v
* output paramtric v
"""
function parametric_divergence_unit_sphere(v)
  function tmp(θϕ)
     Gθϕ=GridapGeosciences.G_unit_sphere(θϕ)
     Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
     function f(θϕ)
        sqrt(det(Gθϕ))*(v(θϕ))
     end
       1.0/sqrt(det(Gθϕ))*(∇⋅f)(θϕ)
  end
end


v(θϕ) = VectorValue(θϕ[1],θϕ[2])
v(θϕ)
parametric_divergence_unit_sphere(v)(θϕ)


divu(θϕ) = parametric_divergence_unit_sphere(_∇u)(θϕ)
divu(θϕ)


"""
* input parametric v
* output paramtric v
"""
function parametric_laplacian_unit_sphere(v)
  function tmp(θϕ)
    parametric_divergence_unit_sphere(parametric_gradient_unit_sphere(v))(θϕ)
  end
end

parametric_laplacian_unit_sphere(_uθϕ)(θϕ)
