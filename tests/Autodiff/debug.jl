using DrWatson
using Gridap
using GridapGeosciences
using Test



my_forward_jacobian(p) = αβ -> transpose( gradient(forward_map(p))(αβ) )

function f(x)
  a = -π/4 + (π/4 - -π/4)*rand()
  b = -π/4 + (π/4 - -π/4)*rand()
  Point(a,b)
end

pts = [f(1) for i in 1:10]
for p in collect(1:6)
  J1 = lazy_map(x->my_forward_jacobian(p)(x), pts)
  J2 = lazy_map(x->forward_jacobian(p)(x), pts)
  @test all(J1 .≈ J2)
end


########### psuedo inverse
function pinvJ(J::Gridap.TensorValues.MultiValue{Tuple{D1,D2}}) where {D1,D2}
  Jt = transpose(J)
  inv(Jt⋅J)⋅Jt
end

my_pinv_jacobian(p) = αβ -> pinvJ(my_forward_jacobian(p)(αβ) )
my_pinv_jacobian(1)(αβ) ≈ forward_pinv_jacobian(1)(αβ)
for p in collect(1:6)
  J1 = lazy_map(x->my_pinv_jacobian(p)(x), pts)
  J2 = lazy_map(x->forward_pinv_jacobian(p)(x), pts)
  @test all(J1 .≈ J2)
end



############# metrics
function metric(p::Int)
  function _m(αβ::Point)
    J = my_forward_jacobian(p)(αβ)
    Jt = transpose(J)
    Jt⋅J
  end
end

for p in collect(1:6)
  m1 = lazy_map(x->metric(p)(x), pts)
  m2 = lazy_map(x->analytic_metric(x), pts)
  @test all(m1 .≈ m2)
end


metric_inv(p) = αβ -> inv(metric(p)(αβ))
for p in collect(1:6)
  m1 = lazy_map(x->metric_inv(p)(x), pts)
  m2 = lazy_map(x->analytic_inv_metric(x), pts)
  @test all(m1 .≈ m2)
end

function perp_matrix(p::Int)
  function _pm(αβ::Point)
    m = metric(p)(αβ)
    TensorValue{2,2}( -m[1,2], m[1,1], -m[2,2], m[1,2] )
  end
end

for p in collect(1:6)
  m1 = lazy_map(x->perp_matrix(p)(x), pts)
  m2 = lazy_map(x->analytic_perp_matrix(x), pts)
  @test all(m1 .≈ m2)
end

my_E(αβ) =  metric(1)(αβ)[1]
my_E(αβ)
E(αβ)

my_detg(p) = αβ -> det(metric(p)(αβ))
my_detg(p)(αβ)

for p in collect(1:6)
  m1 = lazy_map(x->my_detg(p)(x), pts)
  m2 = lazy_map(x->detg(x), pts)
  @test all(m1 .≈ m2)
end

my_sqrtg(p) = αβ -> sqrt(my_detg(p))(αβ)
for p in collect(1:6)
  m1 = lazy_map(x->my_sqrtg(p)(x), pts)
  m2 = lazy_map(x->sqrtg(x), pts)
  @test all(m1 .≈ m2)
end

my_grad_meas(p) = αβ -> gradient(my_sqrtg(p) )(αβ)
for p in collect(1:6)
  m1 = lazy_map(x->my_grad_meas(p)(x), pts)
  m2 = lazy_map(x->grad_meas(x), pts)
  # @test all(m1 .≈ m2)
end
my_grad_meas(p)(αβ)
grad_meas(αβ)
