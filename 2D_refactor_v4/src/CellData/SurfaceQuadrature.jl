
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

function SurfaceQuadrature(m::Metric,quad::CellQuadrature)
  g = m.metric
  SurfaceQuadrature(g,quad)
end

function SurfaceQuadrature(model::ManifoldDiscreteModel,quad::CellQuadrature)
  m = Metric(model)
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

  bgx = lazy_map(Broadcasting(*), bx, gx_meas) # lazy_map(LazyMult(), bx,  gx_meas)

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



struct SurfaceMeasure{C<:SurfaceQuadrature} <: Measure
  s_quad :: C
end

function SurfaceMeasure(s_quad::SurfaceQuadrature)
  C = typeof(s_quad)
  return SurfaceMeasure{C}(s_quad)
end

function SurfaceMeasure(metric,quad::CellQuadrature)
  s_quad = SurfaceQuadrature(metric,quad)
  return SurfaceMeasure(s_quad)
end

function SurfaceMeasure(metric,args...;kwargs...)
  s_quad = SurfaceQuadrature(metric,args...;kwargs...)
  return SurfaceMeasure(s_quad)
end

Gridap.CellData.Measure(s_quad::SurfaceQuadrature) = SurfaceMeasure(s_quad)
# Gridap.CellData.Measure(metric::Function,args...;kwargs...) = SurfaceMeasure(metric,args...;kwargs...)
# Gridap.CellData.Measure(metric::CellField,args...;kwargs...) = SurfaceMeasure(metric,args...;kwargs...)
Gridap.CellData.Measure(metric::Metric,args...;kwargs...) = SurfaceMeasure(metric,args...;kwargs...)

Gridap.CellData.get_cell_quadrature(a::SurfaceMeasure) = a.s_quad.quad

function Gridap.CellData.integrate(f,b::SurfaceMeasure)
  c = integrate(f,b.s_quad)
  cont = DomainContribution()
  add_contribution!(cont,b.s_quad.quad.trian,c)
  cont
end
