using Gridap
using Test
include("../src/initialise.jl")

"""
The theoretical values in _gradf and _divv were computed by hand.
See cubed_sphere_refactor, section 8 notes
"""

function test_operators(m::Metric,
  cf::CellField,_gradf::AbstractArray,
  cv::CellField,_divv::AbstractArray,
  ch::CellField,_laph::AbstractArray,
  pts::CellPoint,xs::AbstractArray
)

  gradf = surface_gradient(cf,m)
  @test sum(gradf(pts) .≈ _gradf) == 1

  divv = surface_divergence(cv,m)
  @test sum(divv(pts) .≈ _divv) == 1

  laph = surface_laplacian(ch,m)
  @test sum((laph(pts)) .≈ _laph) == 1

  # auto diff
  for i in 1:length(xs)
    @test surface_gradient(f,m)(xs[i]) ≈ _gradf[1][i]
    @test surface_divergence(v,m)(xs[i]) ≈ _divv[1][i]
    # @test surface_laplacian(h,m)(xs[i]) ≈ _laph[1][i]
  end

  println("yay!")
end


################################################################################
#### 1D model: x ∈ [0,1]
#### Consider mappings ϕ(x) = (x,f(x))
################################################################################
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1),(1)))
Ω = Triangulation(model)
pts = get_cell_points(Ω)
xs = [Point(0.0,),Point(1.0,)]

f(x) = x[1] # scalar function
v(x) = VectorValue(x[1]) # vector function
h(x) = x[1]^2

cf = CellField(f,Ω)
cv = CellField(v,Ω)
ch = CellField(h,Ω)

### linear mapping, φ(x) = (x,2x)
metric_func(x) = TensorValue{1}(4)
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(1/4,),VectorValue(1/4,)]]
_divv = [[1,1]]
_laph = [[0.5,0.5]]

test_operators(m,cf,_gradf,cv,_divv,ch,_laph,pts,xs)

#### auto diff not implemented
surface_divergence(surface_gradient(h,m),m)(xs[1])


### quadratic mapping, ϕ(x)= (x,x^2)
metric_func(x) = TensorValue{1}(4*x[1]^2 + 1)
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(1,),VectorValue(1/5,)]]
_divv = [[1,9/5]]
_laph = [[2,2/25]]

test_operators(m,cf,_gradf,cv,_divv,ch,_laph,pts,xs)


### cubic mapping, ϕ(x)= (x,x^3 + x^2 + x + 1)
metric_func(x) = TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(1/2,),VectorValue(1/37,)]]
_divv = [[1,85/37]]
_laph = [[1,-22/(37^2)]]

test_operators(m,cf,_gradf,cv,_divv,ch,_laph,pts,xs)

################################################################################
#### 2D model: x,y ∈ [0,1]^2
#### Consider mappings σ(x,y) = (x,y,f(x,y))
################################################################################
model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(1,1)))
Ω = Triangulation(model)
pts = get_cell_points(Ω)
xs = [Point(0.0,0.0),Point(1.0,0.0),Point(0.0,1.0),Point(1.0,1.0)]

f(x) = x[1]^2 + x[1]*x[2] # scalar function
v(x) = VectorValue(x[1]*x[2],x[2]^2) # vector valued function
h(x) = x[1]^2 + x[1]*x[2]^2
cf = CellField(f,Ω)
cv = CellField(v,Ω)
ch = CellField(h,Ω)

### constant mapping, φ(x,y) = (x,y,0)
metric_func(x) = TensorValue{2,2}(1,0,0,1 )
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(0,0),VectorValue(2,1),VectorValue(1,0),VectorValue(3,1)]]
_divv = [[0,0,3,3]]
_laph = [[2,4,2,4]]

test_operators(m,cf,_gradf,cv,_divv,ch,_laph,pts,xs)


### linear mapping, φ(x,y) = (x,y,x+y)
metric_func(x) = TensorValue{2,2}(2,1,1,2 )
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(0,0),VectorValue(1,0),VectorValue(2/3,-1/3),VectorValue(5/3,-1/3)]]
_divv = [[0,0,3,3]]
_laph = [[4/3,8/3,0,4/3]]

test_operators(m,cf,_gradf,cv,_divv,ch,_laph,pts,xs)


### quadratic mapping, φ(x,y) = (x,y,x^2 + y^2)
metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(0,0),VectorValue(2/5,1),VectorValue(1,0),VectorValue(11/9,-7/9)]]
_divv = [[0,0,19/5,35/9]]
_laph = [[2,12/25,2,-164/81]]

test_operators(m,cf,_gradf,cv,_divv,ch,_laph,pts,xs)
