using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Helpers
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools
import Gridap.Helpers: @check
import Gridap.TensorValues: meas

a = 1.0
r = a*sqrt(3.0)

include("maps/metric_maps.jl")
include("surface_metric.jl")
include("surface_operators.jl")

order = 4

parametric_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(8,8))
Ω = Triangulation(parametric_model)
dΩ = Measure(Ω,order)

quad = CellQuadrature(Ω,order)


struct SurfaceQuadrature{DDS,IDS} <: CellDatum
  metric::CellField
  gx::AbstractVector
  gx_meas::AbstractVector
  quad::CellQuadrature{DDS,IDS}
end


function SurfaceQuadrature(g::CellField,quad::CellQuadrature{DDS,IDS}) where {DDS,IDS}

  trian_x = get_triangulation(quad)
  trian_g = get_triangulation(g)

  msg = """\n
    Your are trying to integrate a CellField using a CellQuadrature defined on incompatible
    triangulations. Verify that either the two objects are defined in the same triangulation
    or that the triangulaiton of the CellField is the background triangulation of the CellQuadrature.
    """
  @check is_change_possible(trian_g,trian_x) msg

  # evaluate the metric on the quadtrature domain
  _g = change_domain(g,quad.trian,quad.data_domain_style)
  x = get_cell_points(quad)
  gx = _g(x)
  gx_meas = lazy_map(MetricMeasure(),gx) # get sqrt(meas(g)) for the area element

  SurfaceQuadrature{DDS,IDS}(g,gx,gx_meas,quad)
end


function SurfaceQuadrature(metric::Function,quad::CellQuadrature)
  g = CellField(metric,quad.trian)
  SurfaceQuadrature(g,quad)
end

function SurfaceQuadrature(m::MetricInfo,quad::CellQuadrature)
  g = m.metric
  SurfaceQuadrature(g,quad)
end

function SurfaceQuadrature(metric,args...;kwargs...)
  quad = CellQuadrature(args...;kwargs...)
  SurfaceQuadrature(metric,quad)
end

function Gridap.Fields.integrate(a,s_quad::SurfaceQuadrature)
  println("surface integrate 1")

  quad = s_quad.quad
  b = CellField(a,quad.trian,quad.data_domain_style)
  Gridap.Fields.integrate(b,s_quad)
end

function Gridap.Fields.integrate(f::CellField,s_quad::SurfaceQuadrature)
  println("surface integrate")

  quad = s_quad.quad
  gx_meas = s_quad.gx_meas

  trian_f = get_triangulation(f)
  trian_x = get_triangulation(quad)

  msg = """\n
    Your are trying to integrate a CellField using a CellQuadrature defined on incompatible
    triangulations. Verify that either the two objects are defined in the same triangulation
    or that the triangulaiton of the CellField is the background triangulation of the CellQuadrature.
    """
  @check is_change_possible(trian_f,trian_x) msg

  b = change_domain(f,quad.trian,quad.data_domain_style)
  x = get_cell_points(quad)
  bx = b(x)

  bgx = lazy_map(LazyMult(), bx,  gx_meas)

  if quad.data_domain_style == ReferenceDomain() &&
            quad.integration_domain_style == PhysicalDomain()
    cell_map = get_cell_map(quad.trian)
    cell_Jt = lazy_map(∇,cell_map)
    cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
    lazy_map(IntegrationMap(),bgx,quad.cell_weight,cell_Jtx)
    # lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx,gx)
  else
    @notimplemented
  end
end


# function Gridap.Fields.evaluate!(cache,k::IntegrationMap,aq::AbstractVector,w,jq::AbstractVector,gq::AbstractVector)

#   T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(det(testitem(gq)))
#             + testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(det(testitem(gq))) )
#   z = zero(T)
#   @check length(aq) == length(w)
#   @check length(aq) == length(jq)
#   @inbounds for i in eachindex(aq)
#     z += aq[i]*w[i]*meas(jq[i])*sqrt(det(gq[i]))
#   end
#   z
# end

s_quad = SurfaceQuadrature(metric_func,quad)
out = integrate(1.0,s_quad)
sum(out)*6
4*π*r^2

metric_cf = CellField(metric_func,Ω)
_s_quad = SurfaceQuadrature(metric_cf,quad)
_out = integrate(1.0,_s_quad)
sum(_out)*6

m = MetricInfo(metric_func,Ω)
_s_quad = SurfaceQuadrature(m,quad)
_out = integrate(1.0,_s_quad)
sum(_out)*6


struct SurfaceMeasure{C<:SurfaceQuadrature} <: Measure
  s_quad :: C
end

function SurfaceMeasure(q::SurfaceQuadrature)
  C = typeof(q)
  return SurfaceMeasure{C}(q)
end

function SurfaceMeasure(metric,quad::CellQuadrature)
  s_quad = SurfaceQuadrature(metric,quad)
  return SurfaceMeasure(s_quad)
end

function SurfaceMeasure(metric,args...;kwargs...)
  s_quad = SurfaceQuadrature(metric,args...;kwargs...)
  return SurfaceMeasure(s_quad)
end

Gridap.CellData.Measure(q::SurfaceQuadrature) = SurfaceMeasure(q)
Gridap.CellData.Measure(metric,args...;kwargs...) = SurfaceMeasure(metric,args...;kwargs...)

Gridap.CellData.get_cell_quadrature(a::SurfaceMeasure) = a.s_quad.quad

function Gridap.CellData.integrate(f,b::SurfaceMeasure)
  c = integrate(f,b.s_quad)
  cont = DomainContribution()
  add_contribution!(cont,b.s_quad.quad.trian,c)
  cont
end


dΩg = Measure(metric_func,Ω,order)
out = sum( integrate(1.0,dΩg))
out*6
4*π*r^2

_dΩg = Measure(metric_cf,Ω,order)
_out = sum( integrate(1.0,_dΩg))
_out*6

DΩg = Measure(s_quad)
_out = sum( integrate(1.0,DΩg))
_out*6

dΩg = Measure(m,Ω,order)
out = sum( integrate(1.0,dΩg))
out*6
4*π*r^2
