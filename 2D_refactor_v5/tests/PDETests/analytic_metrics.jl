"""
Consider various mappings:
  1D interval -> polynomial: Ω = [0,1]
    * linear: φ(x) = (x,2x)
    * quad:   ϕ(x) = (x,x^2)
    * cubic:  ϕ(x) = (x,x^3 + x^2 + x + 1)

  2D square -> 3D plane: Ω = [0,1]^2
    * const:  φ(x,y) = (x,y,0)
    * linear: φ(x,y) = (x,y,x+y)
    * quad:   φ(x,y) = (x,y,x^2 + y^2)

  2D square -> cylinder: Ω = [0,2π] × [0,1]
    * 1 periodic boundary in parametric space

  2D square -> sphere: Ω = [0,2π] × [-π/2,π/2]
    * 2 periodic boundaries, requires zeromean constraint and function
    * σ(u,v) = ( cos(u)cos(v), cos(v)sin(u), sin(v) ), u ∈ [0,2π], v ∈ [-π/2,π/2]
"""

using Gridap

################################################################################
#### Metrics
################################################################################

metrics =  Dict{Symbol,Any}()


circle(x) = TensorValue{1}(r^2)
sphere(x) = TensorValue{2,2}(r^2*(cos(x[2]))^2, 0.0, 0.0, r^2 )


linear_1D(x) = TensorValue{1}(4)
quad_1D(x) = TensorValue{1}(4*x[1]^2 + 1)
cubic_1D(x) = TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)

const_2D(x) = TensorValue{2,2}(1.0, 0.0, 0.0, 1.0 )
linear_2D(x) = TensorValue{2,2}(2.0, 1.0, 1.0, 2.0 )
quad_2D(x) = TensorValue{2,2}(1+4*x[1]^2, 4*x[1]*x[2], 4*x[1]*x[2], 1+4*x[2]^2 )




metrics[:linear1D] = linear_1D
metrics[:quad1D] = quad_1D
metrics[:cubic1D] = cubic_1D
metrics[:circle] = circle


metrics[:const2D] = const_2D
metrics[:linear2D] = linear_2D
metrics[:quad2D] = quad_2D
metrics[:cylinder] = const_2D
metrics[:sphere] = sphere

manifolds_1D = [:linear1D,:quad1D,:cubic1D]
manifolds_2D = [:const2D,:linear2D,:quad2D]
manifolds_periodic = [:circle,:sphere]

################################################################################
#### Charts / maps from parametric -> ambient space
################################################################################

map_linear_1D(x) = VectorValue(x[1],2*x[1])
map_quad_1D(x) = VectorValue(x[1],x[1]^2)
map_cubic_1D(x) = VectorValue(x[1],x[1]^3+x[1]^2+x[1]+1)

map_const_2D(x) = VectorValue(x[1],x[2],0.0)
map_linear_2D(x) = VectorValue(x[1],x[2],x[1]+x[2])
map_quad_2D(x) = VectorValue(x[1],x[2],x[1]^2+x[2]^2)

map_circle(x) = VectorValue(r*sin(x[1]),r*cos(x[1]))
map_cylinder(x) = VectorValue(r*sin(x[1]),r*cos(x[1]),x[2])
map_sphere(x) = VectorValue(r*cos(x[1])*cos(x[2]),r*sin(x[1])*cos(x[2]),r*sin(x[2]) )

charts =  Dict{Symbol,Any}()

charts[:linear1D] = map_linear_1D
charts[:quad1D] = map_quad_1D
charts[:cubic1D] = map_cubic_1D
charts[:circle] = map_circle
charts[:const2D] = map_const_2D
charts[:linear2D] = map_linear_2D
charts[:quad2D] = map_quad_2D
charts[:cylinder] = map_const_2D
charts[:sphere] = map_sphere

################################################################################
#### Domains in parametric space
################################################################################

domains =  Dict{Symbol,Any}()
domains[:linear1D] = (0,1)
domains[:quad1D] = (0,1)
domains[:cubic1D] = (0,1)
domains[:circle] = (0,2*π)

domains[:const2D] = (0,1,0,1)
domains[:linear2D] = (0,1,0,1)
domains[:quad2D] = (0,1,0,1)
domains[:cylinder] = (0,2*π,0,1.0)

domains[:sphere] = (0,2*π, -π/2,π/2)


################################################################################
#### Analytic functions periodic on [0,2π] that have
#### ∫ u = ∫ Δu = 0
################################################################################

u1_periodic(x) = cos(x[1])
function u2_periodic(x)
  if x[1] < π
    return x[1]*(π-x[1])
  else
    return (x[1]-π)*(x[1]-2*π)
  end
end

uex_periodic_funcs = Dict{Symbol,Any}()
uex_periodic_funcs[:u1] = u1_periodic
uex_periodic_funcs[:u2] = u2_periodic
