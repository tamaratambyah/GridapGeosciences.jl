using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.TensorValues
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

import Gridap.TensorValues: meas
include("helpers.jl")
include("surface_metric.jl")

a = 1.0
r = a*sqrt(3.0)

order = 4

metric(x) = TensorValue{2,2}(E(x),F(x),F(x),G(x))


parametric_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(8,8))
Ω = Triangulation(parametric_model)
dΩ = Measure(Ω,order)

f = CellField(1.0,Ω)
_metric = CellField(metric,Ω)

quad = CellQuadrature(Ω,order)

b = change_domain(f,quad.trian,quad.data_domain_style)
g = change_domain(_metric,quad.trian,quad.data_domain_style)

x = get_cell_points(quad)

bx = b(x)
gx = g(x)

cell_map = get_cell_map(quad.trian)
cell_Jt = lazy_map(∇,cell_map)
cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)

## compute the integral by hand
weights = collect1d(quad.cell_weight)
jtx = collect1d(cell_Jtx)
sgx = map(x-> sqrt.(meas.(x)), gx)
# _bx = collect1d(lazy_map(SquareMeasure(),TT_x,bx))
z = 0.0

for j in 1:num_cells(parametric_model)
  aq = bx[j]
  jq = jtx[j]
  w = weights[j]
  d = sgx[j]
  @inbounds for i in eachindex(aq)
    z+=(aq[i]*w[i]*(Gridap.TensorValues.meas(jq[i]))*d[i] ) # multiply by sqrt(det(g))
  end

end
z
6*z # assume each face of cube has same area
4*π*r^2

@benchmark map(x-> sqrt.(meas.(x)), gx)

typeof(gx[1]) <: Vector{<:TensorValue}

struct SquareRootMeasure <: Map
end

function Gridap.Arrays.return_cache(f::SquareRootMeasure,cellx::Vector{<:TensorValue})
  x = first(cellx)
  T = typeof((meas(x)))
  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SquareRootMeasure,cellx::Vector{<:TensorValue})
  y = cache

  map!(x->(meas(x)), y, cellx)

  # y .= sqrt.( meas.(cellx)  )

  return y
end


z = lazy_map(SquareRootMeasure(),gx)
cache = array_cache(z)
bm() = lazy_collect(cache,z)
@benchmark bm()



#################################
import Gridap.Helpers: @check
import Gridap.TensorValues: meas

function Gridap.Fields.evaluate!(cache,k::IntegrationMap,aq::AbstractVector,w,jq::AbstractVector,gq::AbstractVector)
  println("integration map")
  T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq)))
            + testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq))) )
  z = zero(T)
  @check length(aq) == length(w)
  @check length(aq) == length(jq)
  @inbounds for i in eachindex(aq)
    z += aq[i]*w[i]*meas(jq[i])#*sqrt(meas(gq[i]))
  end
  z
end

# function Gridap.Fields.return_value(k::IntegrationMap,aq::AbstractArray,w,jq::AbstractVector,gq::AbstractVector)
#   if size(aq,1) == length(w) && size(aq,1) == length(jq)
#     evaluate(k,aq,w,jq,gq)
#   else
#     c = return_cache(k,aq,w,jq,gq)
#     c.array
#   end
# end

# function Gridap.Fields.return_cache(k::IntegrationMap,aq::AbstractArray,w,jq::AbstractVector,gq::AbstractVector)
#   T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq)))
#             + testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq))) )
#   r = zeros(T,size(aq)[2:end])
#   CachedArray(r)
# end

# function Gridap.Fields.evaluate!(cache,k::IntegrationMap,aq::AbstractArray,w,jq::AbstractVector,gq::AbstractVector)
#   setsize!(cache,size(aq)[2:end])
#   r = cache.array
#   @check size(aq,1) == length(w) || size(aq,1) == 0
#   @check size(aq,1) == length(jq) || size(aq,1) == 0
#   fill!(r,zero(eltype(r)))
#   cis = CartesianIndices(r)
#   @inbounds for p in 1:length(w)
#     dV = meas(jq[p])*w[p]*sqrt(meas(gq[p]))
#     for j in cis
#       r[j] += aq[p,j]*dV
#     end
#   end
#   r
# end

# function Gridap.Fields.return_value(k::IntegrationMap,aq::AbstractArray{S,3} where S,w,jq::AbstractVector,gq::AbstractVector)
#   T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq)))
#             + testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq))) )
#   r = zeros(T,size(aq)[2:end])
#   r
# end

# function Gridap.Fields.return_cache(k::IntegrationMap,aq::AbstractArray{S,3} where S,w,jq::AbstractVector,gq::AbstractVector)
#   T = typeof( testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq)))
#             + testitem(aq)*testitem(w)*meas(testitem(jq))*sqrt(meas(testitem(gq))) )
#   r = zeros(T,size(aq)[2:end])
#   s = zeros(typeof(meas(testitem(jq))),length(jq))
#   CachedArray(r), CachedArray(s)
# end

# function Gridap.Fields.evaluate!(cache,k::IntegrationMap,aq::AbstractArray{S,3} where S, w,jq::AbstractVector,gq::AbstractVector)
#   cache_r, cache_s = cache
#   np, ni, nj = size(aq)
#   setsize!(cache_r,(ni,nj))
#   setsize!(cache_s,(np,))
#   r = cache_r.array
#   dV = cache_s.array
#   @check np == length(w) || np == 0
#   @check np == length(jq) || np == 0
#   @inbounds for p in 1:np
#     dV[p] = meas(jq[p])*w[p]*sqrt(meas(gq[p]))
#   end
#   #fill!(r,zero(eltype(r)))
#   @inbounds for j in 1:nj
#     for i in 1:ni
#       rij = zero(eltype(aq))
#       for p in 1:np
#         rij += aq[p,i,j]*dV[p]
#       end
#       r[i,j] = rij
#     end
#   end
#   r
# end

# function Gridap.Fields.evaluate!(cache,k::IntegrationMap,aq::AbstractMatrix, w,jq::AbstractVector,gq::AbstractVector)
#   np, ni = size(aq)
#   setsize!(cache,(ni,))
#   r = cache.array
#   @check np == length(w) || np == 0
#   @check np == length(jq) || np == 0
#   fill!(r,zero(eltype(r)))
#   @inbounds for p in 1:np
#     dV = meas(jq[p])*w[p]*sqrt(meas(gq[p]))
#     for i in 1:ni
#       r[i] += aq[p,i]*dV
#     end
#   end
#   r
# end


out = lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx,gx)
cache = array_cache(out)
bm() = lazy_collect(cache,out)
@benchmark bm()


sum(out)*6
