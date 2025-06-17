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
