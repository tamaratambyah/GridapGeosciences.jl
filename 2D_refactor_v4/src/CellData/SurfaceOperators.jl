function surface_gradient(a::CellField,m::Metric)
  m.inv_metric⋅ gradient(a)
end

function surface_divergence(v::CellField,m::Metric)
  f = m.sq_meas*v
  1/m.sq_meas * divergence(f)
end

function surface_laplacian(f::CellField,m::Metric)
  g = surface_gradient(f,m)
  surface_divergence(g,m)
end


## auto diff
function surface_gradient(f::Number,m::Metric)
  function grad_f(x::Gridap.Fields.Point)
    zero(return_type(outer,x,f))
  end
end


function surface_gradient(f::Function,m::Metric)
  function _gradient(x)
    m.inv_metric(x) ⋅ gradient(f,x)
  end
end

function surface_divergence(f::Function,m::Metric)
  function _divergence(x)
    _f(y) =  m.sq_meas_func(y)*f(y)
    1/m.sq_meas(x) * divergence(_f,x)
  end
end

function surface_laplacian(f::Function,m::Metric)
  @notimplemented
end
