using GridapGeosciences
using Gridap

abstract type GridapGeoscienceType end

struct CubedSphereDescriptor{A<:Real,B<:Real,F<:Function} <: GridapGeoscienceType
  dimension::Int
  radius::A
  thickness::B
  forwardmap::F
end

function CubedSphereDescriptor(dimension::Int,radius=1.0,thickness=1.0)
  forwardmap = forward_map_2D

  if dimension == 3
    forwardmap = forward_map_3D
  end
  A = typeof(radius)
  B = typeof(thickness)
  F = typeof(forwardmap)
  CubedSphereDescriptor{A,B,F}(dimension,radius,thickness,forwardmap)
end

desc = CubedSphereDescriptor(2)
f = desc.forwardmap
f(1)(Point(1,1))

forward_maps

forward_maps = map(p->1,collect(1:6))


function func(p::Int,αβ,radius::Float64)
  @assert length(αβ) == 2
  p*radius
end

function func(radius::Float64)
  function _f(p::Int,αβ)
    func(p,αβ,radius)
  end
end

function func(radius::Float64)
  function _f(p::Int)
    function _g(αβ)
      func(p,αβ,radius)
    end
  end
end

func(1.0)(1)(Point(1,1))
f = func(1.0)
f(1)(Point(1,1))



J(p::Int,x) = jacobians[p](x)
Jt(p::Int,x) =  jacobians_transpose[p](x)
metric(p::Int,x) = Jt(p,x)⋅J(p,x)
inv_metric(p::Int,x) = inv(metric(p,x))
detg(p::Int,x) = det(metric(p,x))
sqrtg(p::Int,x) = sqrt(detg(p,x))
