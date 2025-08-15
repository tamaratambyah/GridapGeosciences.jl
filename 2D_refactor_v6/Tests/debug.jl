αβ = Point(1.0,1.0)
X = forward_map(αβ,1)

_forward_map(p) = αβ -> forward_map(αβ,p)

gradient(_forward_map(1))(αβ)
forward_jacobian(αβ,1)
