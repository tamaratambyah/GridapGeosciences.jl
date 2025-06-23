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

  2D square -> sphere: Ω = [-π/2, π/2 ] × [0,2π]
    * 2 periodic boundaries, requires zeromean constraint and function
"""

using Gridap

metrics =  Dict{Symbol,Any}()


circle(x) = TensorValue{1}(r^2)
sphere(x) = TensorValue{2,2}(r^2, 0.0, 0.0, r^2*(cos(x[1]))^2 )


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

metrics_1D = [:linear1D,:quad1D,:cubic1D]
metrics_2D = [:const2D,:linear2D,:quad2D]


domains =  Dict{Symbol,Any}()
domains[:d1] = (0,1)
domains[:d2] = (0,1,0,1)
domains[:circle] = (0,2*π)
domains[:sphere] = (-π/2,π/2, 0,2*π )
