using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers
using LinearAlgebra

global RADIUS = 1.0

include("forward_map.jl")
include("inverse_map.jl")

p = 1
αβ = Point(π/4,π/4)
αβ_pts = [Point(-π/4,-π/4),Point(0.0,0.0),Point(π/4,π/4)]


forward_map(αβ,p)
forward_jacobian(αβ,p)

f = ForwardMap(p)

XYZ = evaluate(f,αβ)
XYZ_pts = evaluate(f,αβ_pts)

evaluate(∇(f),αβ) == transpose(forward_jacobian(αβ,p))
evaluate(∇(f),pts)


inverse_map(XYZ,p)
inverse_jacobian(XYZ,p)

g = InverseMap(p)

evaluate(g,XYZ)
evaluate(g,XYZ_pts)

evaluate(∇(g),XYZ) == transpose(inverse_jacobian(XYZ,p))
evaluate(∇(g),XYZ_pts)
