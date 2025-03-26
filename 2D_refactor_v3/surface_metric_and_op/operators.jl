function surface_gradient(a::CellField,m::MetricInfo)
  m.inv_metric⋅ gradient(a)
end

function surface_divergence(v::CellField,m::MetricInfo)
  f = m.sq_meas*v
  1/m.sq_meas * divergence(f)
end

function surface_laplacian(f::CellField,m::MetricInfo)
  surface_divergence(surface_gradient(f,m),m)
end


## auto diff
function surface_gradient(f::Number,m::MetricInfo)
  function grad_f(x::Point)
    zero(return_type(outer,x,f))
  end
end


function surface_gradient(f::Function,m::MetricInfo)
  function _gradient(x)
    m.inv_metric(x) ⋅ gradient(f,x)
  end
end

function surface_divergence(f::Function,m::MetricInfo)
  function _divergence(x)
    _f(y) =  m.sq_meas_func(y)*f(y)
    1/m.sq_meas(x) * divergence(_f,x)
  end
end
