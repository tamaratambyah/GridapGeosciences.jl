using Gridap
using Test
include("../src/initialise.jl")

function test_operators(m::Metric,
  cf::CellField,_gradf::AbstractArray,
  cv::CellField,_divv::AbstractArray,
  pts::CellPoint,xs::AbstractArray
)

  gradf = surface_gradient(cf,m)
  @test sum(gradf(pts) .≈ _gradf) == 1
  # println(gradf(pts))

  divv = surface_divergence(cv,m)
  @test sum(divv(pts) .≈ _divv) == 1
  # println(divv(pts))

  # auto diff
  for i in 1:length(xs)
    @test surface_gradient(f,m)(xs[i]) ≈ _gradf[1][i]
    @test surface_divergence(v,m)(xs[i]) ≈ _divv[1][i]
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
cf = CellField(f,Ω)
cv = CellField(v,Ω)


### linear mapping, φ(x) = (x,2x)
metric_func(x) = TensorValue{1}(4)
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(1/4,),VectorValue(1/4,)]]
_divv = [[1,1]]

test_operators(m,cf,_gradf,cv,_divv,pts,xs)



gradf = surface_gradient(cf,m)
@test sum(gradf(pts) .≈ _gradf) == 1
println(gradf(pts))

divv = surface_divergence(cv,m)
@test sum(divv(pts) .≈ _divv) == 1
println(divv(pts))

# auto diff
for i in 1:length(xs)
  @test surface_gradient(f,m)(xs[i]) == _gradf[1][i]
  @test surface_divergence(v,m)(xs[i]) == _divv[1][i]
end

### quadratic mapping, ϕ(x)= (x,x^2)
metric_func(x) = TensorValue{1}(4*x[1]^2 + 1)
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(1,),VectorValue(1/5,)]]
_divv = [[1,9/5]]

gradf = surface_gradient(cf,m)
@test sum(gradf(pts) .≈ _gradf) == 1
println(gradf(pts))

divv = surface_divergence(cv,m)
@test sum(divv(pts) .≈ _divv) == 1
println(divv(pts))

# auto diff
for i in 1:length(xs)
  @test surface_gradient(f,m)(xs[i]) ≈ _gradf[1][i]
  @test surface_divergence(v,m)(xs[i]) ≈ _divv[1][i]
end


### cubic mapping, ϕ(x)= (x,x^3 + x^2 + x + 1)
metric_func(x) = TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(1/2,),VectorValue(1/37,)]]
_divv = [[1,85/37]]

gradf = surface_gradient(cf,m)
@test sum(gradf(pts) .≈ _gradf) == 1
println(gradf(pts))

divv = surface_divergence(cv,m)
@test sum(divv(pts) .≈ _divv) == 1
println(divv(pts))

# auto diff
for i in 1:length(xs)
  @test surface_gradient(f,m)(xs[i]) ≈ _gradf[1][i]
  @test surface_divergence(v,m)(xs[i]) ≈ _divv[1][i]
end

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
cf = CellField(f,Ω)
cv = CellField(v,Ω)

### constant mapping, φ(x,y) = (x,y,0)
metric_func(x) = TensorValue{2,2}(1,0,0,1 )
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(0,0),VectorValue(2,1),VectorValue(1,0),VectorValue(3,1)]]
_divv = [[0,0,3,3]]

gradf = surface_gradient(cf,m)
@test sum(gradf(pts) .≈ _gradf) == 1
println(gradf(pts))

divv = surface_divergence(cv,m)
@test sum(divv(pts) .≈ _divv) == 1
println(divv(pts))

# auto diff
for i in 1:length(xs)
  @test surface_gradient(f,m)(xs[i]) ≈ _gradf[1][i]
  @test surface_divergence(v,m)(xs[i]) ≈ _divv[1][i]
end

### linear mapping, φ(x,y) = (x,y,x+y)
metric_func(x) = TensorValue{2,2}(2,1,1,2 )
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(0,0),VectorValue(1,0),VectorValue(2/3,-1/3),VectorValue(5/3,-1/3)]]
_divv = [[0,0,3,3]]

gradf = surface_gradient(cf,m)
@test sum(gradf(pts) .≈ _gradf) == 1
println(gradf(pts))

divv = surface_divergence(cv,m)
@test sum(divv(pts) .≈ _divv) == 1
println(divv(pts))

# auto diff
for i in 1:length(xs)
  @test surface_gradient(f,m)(xs[i]) ≈ _gradf[1][i]
  @test surface_divergence(v,m)(xs[i]) ≈ _divv[1][i]
end

### quadratic mapping, φ(x,y) = (x,y,x^2 + y^2)
metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )
m = Metric(metric_func,Ω)
_gradf = [[VectorValue(0,0),VectorValue(2/5,1),VectorValue(1,0),VectorValue(11/9,-7/9)]]
_divv = [[0,0,19/5,35/9]]

gradf = surface_gradient(cf,m)
@test sum(gradf(pts) .≈ _gradf) == 1
println(gradf(pts))

divv = surface_divergence(cv,m)
@test sum(divv(pts) .≈ _divv) == 1
println(divv(pts))

# auto diff
for i in 1:length(xs)
  @test surface_gradient(f,m)(xs[i]) ≈ _gradf[1][i]
  @test surface_divergence(v,m)(xs[i]) ≈ _divv[1][i]
end
