"""
Consider surface metric of the form:
g = [E F
     F G]
"""

function factor(x)
  α,β = x
   r^2/( (1 + (tan(α))^2 + (tan(β))^2 )^2 * (cos(α))^2 * (cos(β))^2 )
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

metric_func(x) = TensorValue{2,2}(E(x),F(x),F(x),G(x))
