using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")
coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
model = Adaptivity.refine(coarse_model)
ref_model = Adaptivity.refine(model)

manifold_model = coarse_model
ambient_model = AmbientDiscreteModel(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
get_background_model(Ω_parametric) == get_background_model(Ω_ambient)

####### parametric -> ambient
cf = CellField(x->x[1],Ω_parametric)
V = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1); conformity=:H1)
uh = interpolate(cf,V)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf, "uh"=>uh],append=false)

_cf = change_domain(cf,Ω_ambient,DomainStyle(cf))
_uh = change_domain(uh,Ω_ambient,DomainStyle(uh))
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>_cf, "uh"=>_uh],append=false)


####### ambient -> parametric
cf = CellField(x->x[1],Ω_ambient)
V = FESpace(ambient_model,ReferenceFE(lagrangian,Float64,1); conformity=:H1)
uh = interpolate(cf,V)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf, "uh"=>uh],append=false)

_cf = change_domain(cf,Ω_parametric,DomainStyle(cf))
_uh = change_domain(uh,Ω_parametric,DomainStyle(uh))
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>_cf, "uh"=>_uh],append=false)


###### periodic functions on sphere
# f(x) = cos(2*pi*x[1])*cos(2*pi*x[2])*cos(2*pi*x[3])
f(x) = x[1]*x[2]#*x[3]
cf = CellField(f,Ω_ambient)
V = FESpace(ambient_model,ReferenceFE(lagrangian,Float64,1); conformity=:H1)
uh = interpolate(cf,V)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf, "uh"=>uh],append=false)

ambient_coords = get_cell_coordinates(Ω_ambient)
map(x->f(x),ambient_coords[1])

z = lazy_map(ConvertMap(f),ambient_coords)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)


V_parametric = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1); conformity=:H1)
uh = CellField(V_parametric,collect(z))
writevtk(Ω_parametric,dir*"/parametric",cellfields=["uh"=>uh],append=false)

struct ConvertMap{A} <: Map
  f::A
end

function Gridap.Arrays.return_cache(k::ConvertMap,cellx::AbstractArray)
  x = first(cellx)
  T = typeof(k.f(x))
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,k::ConvertMap,cellx::AbstractArray)
  y = cache
  map!(x->k.f(x), y, cellx)
  return y
end

_cf = change_domain(cf,Ω_parametric,DomainStyle(cf))
_uh = change_domain(uh,Ω_parametric,DomainStyle(uh))
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>_cf, "uh"=>_uh],append=false)
