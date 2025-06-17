using Gridap
using Plots, LaTeXStrings
using DrWatson

include("../../src/initialise.jl")


function wave_divergence(v::CellField,m::Metric)
  println("cellfield ")
  (
   m.sq_meas*(divergence(m.inv_metric)⋅v
            + Operation(tr)( ( m.inv_metric⋅ gradient(v) )) )
  + gradient(m.sq_meas) ⋅ (m.inv_metric ⋅ v)
  )
end

function wave_divergence(v::Function,m::Metric)
  println("function")
  function _wave_divergence(x)
    (
      m.sq_meas_func(x)*( divergence(m.inv_metric_func,x)⋅v(x)
                        + tr( m.inv_metric_func(x) ⋅  gradient(v,x) )
                        )
     + gradient(m.sq_meas_func,x) ⋅ (m.inv_metric_func(x) ⋅ v(x))
     )
  end
end


wave_divergence(f::Function,x::SVector) = wave_divergence(f,Point(x))

metric_func(x) = TensorValue{1}( 1 + 4*x[1]^2 )

model = CartesianDiscreteModel((0,1), (8,),isperiodic=(true,))
Ω = Triangulation(model)
m = Metric(metric_func,Ω)

v(x) = VectorValue(x[1])
pex(x) = x[1]

vcf = CellField(v,Ω)
pcf = CellField(pex,Ω)

pt = Point(1)

pex(pt)*(wave_divergence(v,m))(pt)
(pcf*wave_divergence(vcf,m))(pt)

u(x) = x[1]
_τ(x) = 1.0
τ(x) = VectorValue(_τ(x))
f(x) = (4*x[1]^2+1)^(-0.5)*_τ(x)

uf(x) = -1.0*u(x)*gradient(f)(x)
uf(pt)

_uf(x) = -1.0*u(x)*(wave_divergence(τ,m))(x)
_uf(pt)
