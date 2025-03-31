using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")


n = 10
cell_pts = [Point(1.0,-1.0,-1.0), Point(1.0,1.0,-1.0), Point(1.0,-1.0,1.0), Point(1.0,1.0,1.0)]
cellx = fill(cell_pts,n)
panel_ids = [rand(1:6) for i in 1:n]


################################################################################
### Rotation
################################################################################

### no panel ids
panel = 1
f = PanelRotationField(rp1[panel])
fn = fill(f,n)

evaluate(f,cell_pts[1])
evaluate(f,cell_pts)

z = lazy_map(evaluate,fn,cellx)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)
getindex!(cache,z,1)

### gradients
gradf = gradient(f)

evaluate(gradf,cell_pts[1])

z = lazy_map(gradf,cell_pts)
cache = array_cache(z)
print_lazy_arr(z)
@benchmark lazy_collect(cache,z)


gradfn = lazy_map(Broadcasting(∇), fn)

z = lazy_map(gradfn,cellx)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)
getindex!(cache,z,1)


### with panel ids
f = PanelRotationField(rp1)
fn = map(x->PanelRotationField(rp1[x]),panel_ids)
evaluate(f,cell_pts[1],panel_ids[panel])
evaluate(f,cell_pts,panel_ids[panel])


z = lazy_map(evaluate,fn,cellx,panel_ids)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)

### gradeints
gradf = gradient(f)
evaluate(gradf,cell_pts[1],panel_ids[1])

z = lazy_map(gradf,cell_pts,panel_ids)
cache = array_cache(z)
print_lazy_arr(z)
@benchmark lazy_collect(cache,z)

gradfn = lazy_map(Broadcasting(∇), fn)

z = lazy_map(gradfn,cellx)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)
getindex!(cache,z,1)

#### map tests:
z = lazy_map(PanelMap(),cell_pts,panel_ids)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)


################################################################################
### Bump
################################################################################
g = BumpField(A_bump,B_bump,b_bump)
gn = fill(g,n)

evaluate(g,cell_pts[1])
evaluate(g,cell_pts)

z = lazy_map(evaluate,gn,cellx)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)

## gradients
gradg = gradient(g)
evaluate(gradg,cell_pts[1])

z = lazy_map(gradg,cell_pts)
cache = array_cache(z)
print_lazy_arr(z)
@benchmark lazy_collect(cache,z)

gradgn = lazy_map(Broadcasting(∇), gn)

z = lazy_map(gradgn,cellx)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)
getindex!(cache,z,1)

#### map tests:
z = lazy_map(BumpMap(),cell_pts)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)



################################################################################
#### Composition of fields
################################################################################
f = PanelRotationField(rp1[1])
g = BumpField(A_bump,B_bump,b_bump)

k = g∘f
evaluate(k,cell_pts)

z = lazy_map(k,cellx)
getindex!(array_cache(z),z,1)

gradk = gradient(k)
z = lazy_map(gradk,cellx)
getindex!(array_cache(z),z,1)


# fn = lazy_map(x-> PanelRotationField(rp1[x]), panel_ids)
# gn = fill(g,n)
# kn = lazy_map((∘),gn, fn)
kn = lazy_map(x-> g∘PanelRotationField(rp1[x]), panel_ids  )
z = lazy_map(evaluate,kn,cellx)
cache = array_cache(z)
getindex!(cache,z,1)
@benchmark lazy_collect(cache,z)


gradkn = lazy_map(Broadcasting(∇),kn)
z = lazy_map(evaluate,gradkn,cellx)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)
getindex!(cache,z,1)


################################################################################
#### Composition of Maps
################################################################################

_cell_pts = lazy_map(PanelMap(),cell_pts,panel_ids)
z = lazy_map(BumpMap(),_cell_pts)
cache = array_cache(z)
@benchmark lazy_collect(cache,z)
