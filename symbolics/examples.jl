#inline, errors in REPL
using DrWatson
using Symbolics

@variables x y
@variables t

D = Differential(t)

z = t + t^2
expand_derivatives(D(z))
Symbolics.derivative(z,t;simplify=true)

Symbolics.jacobian([x + x*y, x^2 + y], [x, y])


function f(u)
  x,y = u
  [x + x*y, x^2 + y]
end
u = Symbolics.variables(:u, 1:2)
f(u)
Symbolics.derivative(f(u),u;simplify=true)
Symbolics.jacobian(f(u), u;simplify=true)


function g(u)
  2.0 .*u
end

g(f(u))

Symbolics.jacobian(g(f(u)), u;simplify=true)
