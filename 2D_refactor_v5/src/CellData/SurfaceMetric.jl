"""
Metric

A structure to hold the metric information required for surface operations and
surface quadrature.
The idea is that the user gives a triangulation and metric_function of the form

    metric_func(x) = TensorValue{2,2}(E(x),F(x),F(x),G(x))

Then, we compute the inverse and sqrt(det) as cell fields on the triangulation,
and pass this information to the surface quadrature/operators.

Using metric_func(x), we generate sq_meas_func(x) and inv_metric_func(x). This
assists with automatic differentiation in surface operators.
In addition, we return CellFields for sq_meas and inv_metric which allows us to
mulitply with other cell fields via OperationCellFields.

Note, I tried making sq_meas and inv_metric OperationCellFields from the start,
but this caused issues for auto diff divergence. This is why I created the
functions.
"""

struct Metric{A<:CellField,B<:CellField,C<:CellField} <: CellDatum
  metric::A
  sq_meas::B
  inv_metric::C
  metric_func::Function
  sq_meas_func::Function
  inv_metric_func::Function
end

function Metric(metric::CellField,sq_meas::CellField,inv_metric::CellField,
  metric_func,sq_meas_func,inv_metric_func)
  A = typeof(metric)
  B = typeof(sq_meas)
  C = typeof(inv_metric)
  Metric{A,B,C}(metric,sq_meas,inv_metric,metric_func,sq_meas_func,inv_metric_func)
end

function Metric(metric_func::Function,Ω::Triangulation)
  sq_meas_func(x) = sqrt(meas(metric_func(x)))
  inv_metric_func(x) = inv(metric_func(x))

  metric = CellField(metric_func,Ω)
  sq_meas = CellField(sq_meas_func,Ω)
  inv_metric = CellField(inv_metric_func,Ω)
  Metric(metric,sq_meas,inv_metric,metric_func,sq_meas_func,inv_metric_func)
end

function Metric(::CubedSphere,Ω::Triangulation)
  println("cubed sphere metric")
  _metric_func(x) = metric_func(cubedsphere)(x)
  _sq_meas_func(x) = sq_meas_func(cubedsphere)(x)
  _inv_metric_func(x) = inv_metric_func(cubedsphere)(x)

  metric = CellField(_metric_func,Ω)
  sq_meas = CellField(_sq_meas_func,Ω)
  inv_metric = CellField(_inv_metric_func,Ω)
  Metric(metric,sq_meas,inv_metric,_metric_func,_sq_meas_func,_inv_metric_func)


end


function Metric(name::ManifoldName,Ω::Triangulation)
  _metric_func(x) = metric_func(name)(x)
  Metric(_metric_func,Ω)
end


function Metric(model::ManifoldDiscreteModel)
  manifold_name = get_manifold_name(model)

  _metric_func(x) = metric_func(manifold_name)(x)
  Ω = Triangulation(model)

  Metric(_metric_func,Ω)
end

function Metric(model::AdaptedDiscreteModel)
  manifold_name = get_manifold_name(model)

  # _metric_func(x) = metric_func(manifold_name)(x)
  Ω = Triangulation(model)


  println("adapted cubed sphere metric")
  _metric_func(x) = metric_func(cubedsphere)(x)
  _sq_meas_func(x) = sq_meas_func(cubedsphere)(x)
  _inv_metric_func(x) = inv_metric_func(cubedsphere)(x)

  metric = CellField(_metric_func,Ω)
  sq_meas = CellField(_sq_meas_func,Ω)
  inv_metric = CellField(_inv_metric_func,Ω)


  Metric(metric,sq_meas,inv_metric,_metric_func,_sq_meas_func,_inv_metric_func)
end

"""
Consider surface metric of the form:
g = [E F
     F G]
"""
function metric_func(::CubedSphere)
  function _metric_func(x)
    TensorValue{2,2}(E(x),F(x),F(x),G(x))
  end
end

function inv_metric_func(::CubedSphere)
  function _inv_metric_func(x)
    1/(E(x)*G(x) - F(x)*F(x))*TensorValue{2,2}(G(x),-1.0*F(x),-1.0*F(x),E(x))
  end
end

function det(::CubedSphere)
  function _det(x)
    E(x)*G(x) - F(x)*F(x)
  end
end

function sq_meas_func(::CubedSphere)
  function _sq_meas_func(x)
    sqrt( E(x)*G(x) - F(x)*F(x) )
  end
end

function metric_func(::Cube)
  function _metric_func(x)
    TensorValue{2,2}(1.0,0.0,0.0,1.0)
  end
end

function factor(x)
  α,β = x
  RADIUS^2/( (1 + (tan(α))^2 + (tan(β))^2 )^2 * (cos(α))^2 * (cos(β))^2 )
end



function E(x)
  α,β = x
  factor(x)*( 1 + (tan(α))^2 )
end

function F(x)
  α,β = x
  -1.0*factor(x)*( tan(α)*tan(β)  )
end

function G(x)
  α,β = x
  factor(x)*( 1 + (tan(β))^2 )
end
