"""
Test the new wave divergence operator that arises from the following IBP on
doubely periodic domains:
  ∫( surf_grad(u)⋅τ  )dΩ_g = - ∫( u wave_div(√g G^{-1} τ  ) )dΩ

where
  - u is scalar, τ is vector
  - surf_grad, dΩ_g include the metric
  - wave_div is the new operator
  - dΩ is regular quadrature
  - Ω is doubly periodic

The LHS can be written as:
  ∫( (G^{-1}⋅∇u)⋅τ √g  )dΩ  = ∫( (√g G^{-1}⋅τ) ⋅ ∇u  )dΩ
                            = -∫( u ∇⋅(√g G^{-1}⋅τ ) )dΩ (after IBP; periodic BC)

This is to test that
   ∇⋅(√g G^{-1}⋅τ ) = wave_div(√g G^{-1}τ)

Consider various mappings:
  1D interval -> polynomial: Ω = [0,1]
    * linear: φ(x) = (x,2x)
    * quad:   ϕ(x) = (x,x^2)
    * cubic:  ϕ(x) = (x,x^3 + x^2 + x + 1)

  2D square -> 3D plane: Ω = [0,1]^2
    * const:  φ(x,y) = (x,y,0)
    * linear: φ(x,y) = (x,y,x+y)
    * quad:   φ(x,y) = (x,y,x^2 + y^2)
"""

using Gridap
using Plots, LaTeXStrings
using DrWatson

include("../../src/initialise.jl")

################################################################################
#### Divergence definitions on the LHS and RHS
################################################################################
LHSinside(x) = sqrtg(x)*(Ginv(x)⋅τ(x))
divLHS(x) = divergence(LHSinside)(x)

divRHS(x) = sqrtg(x)*( divergence(Ginv)(x)⋅τ(x) + tr( Ginv(x)⋅ gradient(τ)(x)  )) + gradient(sqrtg)(x)⋅(Ginv(x)⋅τ(x))


################################################################################
#### 1D tests
################################################################################
model = CartesianDiscreteModel((0,1), (8,),isperiodic=(true,))
Ω = Triangulation(model)

pt = Point(1)


τ(x) = VectorValue(x[1]^2+ x[1]^3 + 4)
τcf = CellField(τ,Ω)

### linear mapping, φ(x) = (x,2x)
sqrtg(x) = sqrt( 4  )
Ginv(x) = TensorValue{1}(1/(4) )
@test divLHS(pt) ≈ divRHS(pt)

metric_func(x) = TensorValue{1}(4)
m = Metric(metric_func,Ω)
@test divRHS(pt) ≈  wave_divergence(τ,m)(pt)
@test divRHS(pt) ≈  wave_divergence(τcf,m)(pt)


### quadratic mapping, ϕ(x)= (x,x^2)
sqrtg(x) = sqrt( 1 + 4*x[1]^2  )
Ginv(x) = TensorValue{1}(1/(1+4*x[1]^2) )
@test divLHS(pt) ≈ divRHS(pt)

metric_func(x) = TensorValue{1}(4*x[1]^2 + 1)
m = Metric(metric_func,Ω)
@test divRHS(pt) ≈  wave_divergence(τ,m)(pt)
@test divRHS(pt) ≈  wave_divergence(τcf,m)(pt)


### cubic mapping, ϕ(x)= (x,x^3 + x^2 + x + 1)
sqrtg(x) = sqrt( 9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2  )
Ginv(x) = TensorValue{1}(1/(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2) )
@test divLHS(pt) ≈ divRHS(pt)

metric_func(x) = TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)
m = Metric(metric_func,Ω)
@test divRHS(pt) ≈  wave_divergence(τ,m)(pt)
@test divRHS(pt) ≈  wave_divergence(τcf,m)(pt)


################################################################################
#### 2D tests
################################################################################
model = CartesianDiscreteModel((0,1,0,1), (8,8),isperiodic=(true,true))
Ω = Triangulation(model)

pt = Point(1,1)

τ(x) = VectorValue(x[1]^2*x[2],x[2] + x[1]^3*x[2]^3)
τcf = CellField(τ,Ω)

### constant mapping, φ(x,y) = (x,y,0)
sqrtg(x) = 1.0
Ginv(x) = TensorValue{2,2}( 1.0, 0.0,0.0,1.0   )
@test divLHS(pt) ≈ divRHS(pt)

metric_func(x) = TensorValue{2,2}(1,0,0,1 )
m = Metric(metric_func,Ω)
@test divRHS(pt) ≈  wave_divergence(τ,m)(pt)
@test divRHS(pt) ≈  wave_divergence(τcf,m)(pt)


### linear mapping, φ(x,y) = (x,y,x+y)
sqrtg(x) = sqrt(3)
Ginv(x) = TensorValue{2,2}( 2/3, -1/3, -1/3, 2/3  )
@test divLHS(pt) ≈ divRHS(pt)

metric_func(x) = TensorValue{2,2}(2,1,1,2 )
m = Metric(metric_func,Ω)
@test divRHS(pt) ≈  wave_divergence(τ,m)(pt)
@test divRHS(pt) ≈  wave_divergence(τcf,m)(pt)


### quadratic mapping, φ(x,y) = (x,y,x^2 + y^2)
sqrtg(x) = sqrt( 1 + 4*x[1]^2 + 4*x[2]^2 )
A(x) =  (1+4*x[2]^2)/( 1 + 4*x[1]^2 + 4*x[2]^2 )
B(x) = -4*x[1]*x[2]/( 1 + 4*x[1]^2 + 4*x[2]^2 )
D(x) = (1+4*x[1]^2)/( 1 + 4*x[1]^2 + 4*x[2]^2 )
Ginv(x) = TensorValue{2,2}( A(x),B(x),B(x),D(x)   )

@test divLHS(pt) ≈ divRHS(pt)

metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )
m = Metric(metric_func,Ω)
@test divRHS(pt) ≈  wave_divergence(τ,m)(pt)
@test divRHS(pt) ≈  wave_divergence(τcf,m)(pt)
