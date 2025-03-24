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

include("surface_metric.jl")

a = 1.0
r = a*sqrt(3.0)

order = 4

parametric_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(8,8))
Ω = Triangulation(parametric_model)
dΩ = Measure(Ω,order)

f = CellField(1.0,Ω)
quad = CellQuadrature(Ω,order)

metric(x) = TensorValue{2,2}(E(x),F(x),F(x),G(x))


struct SurfaceQuadrature{DDS,IDS} <: CellDatum
  quad::CellQuadrature{DDS,IDS}
  metric::CellField
end

function SurfaceQuadrature(quad::CellQuadrature{DDS,IDS},metric::Function) where {DDS,IDS}
  g = CellField(metric,quad.trian)
  SurfaceQuadrature{DDS,IDS}(quad,g)
end


s_quad = SurfaceQuadrature(quad,metric)


out = integrate(1,s_quad)
sum(out)*6



import Gridap.Helpers: @check
import Gridap.TensorValues: meas

function Gridap.Fields.integrate(a,s_quad::SurfaceQuadrature)
  println("surface integrate 1")

  quad = s_quad.quad
  b = CellField(a,quad.trian,quad.data_domain_style)
  Gridap.Fields.integrate(b,s_quad)
end

function Gridap.Fields.integrate(f::CellField,s_quad::SurfaceQuadrature)
  println("surface integrate")
  quad = s_quad.quad
  metric = s_quad.metric

  trian_f = get_triangulation(f)
  trian_x = get_triangulation(quad)
  trian_g = get_triangulation(metric)

  msg = """\n
    Your are trying to integrate a CellField using a CellQuadrature defined on incompatible
    triangulations. Verify that either the two objects are defined in the same triangulation
    or that the triangulaiton of the CellField is the background triangulation of the CellQuadrature.
    """
  @check is_change_possible(trian_f,trian_x) msg

  b = change_domain(f,quad.trian,quad.data_domain_style)
  g = change_domain(metric,quad.trian,quad.data_domain_style)
  x = get_cell_points(quad)
  bx = b(x)
  gx = g(x)

  # bgx = lazy_map(SquareMeasure(),gx,bx)

  if quad.data_domain_style == ReferenceDomain() &&
            quad.integration_domain_style == PhysicalDomain()
    cell_map = get_cell_map(quad.trian)
    cell_Jt = lazy_map(∇,cell_map)
    cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
    lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx,gx)
  else
    @notimplemented
  end
end

import Gridap.Helpers: @check

function Gridap.Fields.evaluate!(cache,k::IntegrationMap,aq::AbstractVector,w,jq::AbstractVector,gq::AbstractVector)
  println("integration map")
  T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq)))
            + testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq))) )
  z = zero(T)
  @check length(aq) == length(w)
  @check length(aq) == length(jq)
  @inbounds for i in eachindex(aq)
    z += aq[i]*w[i]*meas(jq[i])*sqrt(meas(gq[i]))
  end
  z
end
