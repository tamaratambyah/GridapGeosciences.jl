"""
MetricInfo

A structure to hold the metric information required for surface operations and
surface quadrature.
The idea is that the user gives a triangulation and metric_function of the form

    metric_func(x) = TensorValue{2,2}(E(x),F(x),F(x),G(x))

Then, we compute the inverse and sqrt(det) as cell fields on the triangulation,
and pass this information to the quadrature/operators.

Use functions to compute inverse and sqrt(det) such that we can return CellFields.
This allows us to mulitply with other cell fields via OperationCellFields.

Note, I tried making sq_meas and inv_metric OperationCellFields from the beginning,
but this caused issues for auto diff divergence. This is why I created the
functions.
"""

struct MetricInfo{A<:CellField,B<:CellField,C<:CellField} <: CellDatum
  metric::A
  sq_meas::B
  inv_metric::C
  metric_func::Function
  sq_meas_func::Function
  inv_metric_func::Function
end

function MetricInfo(metric::CellField,sq_meas::CellField,inv_metric::CellField,
  metric_func,sq_meas_func,inv_metric_func)
  A = typeof(metric)
  B = typeof(sq_meas)
  C = typeof(inv_metric)
  MetricInfo{A,B,C}(metric,sq_meas,inv_metric,metric_func,sq_meas_func,inv_metric_func)
end

function MetricInfo(metric_func::Function,Ω::Triangulation)
  sq_meas_func(x) = sqrt(meas(metric_func(x)))
  inv_metric_func(x) = inv(metric_func(x))

  metric = CellField(metric_func,Ω)
  sq_meas = CellField(sq_meas_func,Ω)
  inv_metric = CellField(inv_metric_func,Ω)
  MetricInfo(metric,sq_meas,inv_metric,metric_func,sq_meas_func,inv_metric_func)
end
