using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Helpers
using Gridap.TensorValues
using Gridap.Fields
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

include("maps/map_matrices.jl")
rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()
rp1 = map(TensorValue,rotate_panel_p_to_1)
r1p = map(TensorValue,rotate_panel_1_to_p)

_A , _B, _b = bump_matrics(1.0)
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)


################################################################################
#### Rotation maps
struct PanelRotationField{A} <: Gridap.Fields.Field
  mats::A
end


function Gridap.Arrays.return_cache(f::PanelRotationField,x::VectorValue)
  y = zero(x)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationField,x::VectorValue)
  println("evaluate f at point")
  y = cache
  A = f.mats
  y = A .⋅ x
  return y
end

mats = rp1[1]
pt = Point(1.0,1.0,1.0)
f = PanelRotationField(mats)

evaluate(f,pt)

## gradient
# the gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:PanelRotationField},x::VectorValue)
  f.object.mats
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:PanelRotationField},x::VectorValue)
  println("evaluate grad f at point")

  w = cache
  w = f.object.mats
  return transpose(w)
end


gradf = gradient(f)
gradf(pt)


################################################################################
#### Bump 3D <-> 2D maps
struct Panel1BumpField{A,B,b} <: Gridap.Fields.Field
  A_bump::A
  B_bump::B
  b_bump::b
end


function Gridap.Arrays.return_cache(f::Panel1BumpField,x::VectorValue{D}) where {D}
  A = f.A_bump
  B = f.B_bump

  if D == 3 # D==3, -> bump 3D -> 2D, y == 2 components;
    T = typeof(A⋅x)
  elseif D == 2 # D==2 -> bump 2D -> 3D, y == 3 components;
    T = typeof(B⋅x)
  end

  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Panel1BumpField,x::VectorValue{D}) where {D}
  println("eval")
  y = cache
  A = f.A_bump
  B = f.B_bump
  b = f.b_bump

  if D == 3 # bump 3D -> 2D
    y = A.⋅x
  elseif D == 2 # bump 2D -> 3D
    y = B.⋅x .+ b
  end

  return y
end

g = Panel1BumpField(A_bump,B_bump,b_bump)

evaluate(g,pt)

k = g∘f
evaluate(k,pt)


### gradient
# the gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:Panel1BumpField},x::VectorValue{D}) where {D}

  if D == 3
    T = typeof(f.object.A_bump)
  elseif D == 2
    T = typeof(f.object.B_bump)
  end

  y = zero(T)

  return y
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:Panel1BumpField},x::VectorValue{D}) where {D}
  println("grad eval")
  y = cache
  println(D)
  if D == 3
   y = f.object.A_bump
  elseif D == 2
    y = f.object.B_bump
  end

  return transpose(y)

end

# # f ∘ g = f(g(x))
# function Base.:∘(f::Field,g::Field)
#   print("my composition")
#   Operation(g)(f)
# end

gradient(g)(pt)

#### composition

# function Gridap.Fields.gradient(f::Fields.OperationField{<:Field})
#   a = f.op
#   @notimplementedif length(f.fields) != 1
#   b, = f.fields
#   x = ∇(a)∘b
#   y = ∇(b)
#   # y⋅x
#   x⋅y
# end

k = g ∘ f
w = gradient(k)
gradk = w(pt)
get_array(gradk)

_k = f ∘ g
_w = gradient(_k)
_gradk = _w(Point(1.0,1.0))
get_array(_w(Point(1.0,1.0)))



# x = Point(1.0,1.0)
# fa(x) = VectorValue(x[1],2.0)
# fb(x) = x[1] + x[2]

# f = GenericField(fa)
# g = GenericField(fb)

# k = g∘f
# k(x)

# gradk = gradient(k)
# gradk(x)

# _k = f∘g
# _k(x)

# _gradk = gradient(_k)
# _gradk(x)
# _fba = Operation(fb)(fa)
# evaluate(_fba,x)

# fab = Operation(fa)(fb)
# c = evaluate(fab,x)


########### compose with cmaps
include("geometry/cube_surface_1_cell_per_panel.jl")
cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel(1.0)...)
Ω = Triangulation(cube_model_3D)
dΩ = Measure(Ω,3)
grid = get_grid(cube_model_3D)

cmap = get_cell_map(grid)
evaluate(cmap[1],Point(1,1))
f
_f = fill(f,num_cells(grid))
_g = fill(g,num_cells(grid))
_k = fill(g ∘ f,num_cells(grid))
_cmap = lazy_map(∘,_k,cmap)
evaluate(_cmap[1],Point(1,1))


quad = dΩ.quad
cell_Jt = lazy_map(∇,_cmap)
cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
map(x->Gridap.TensorValues.meas(x),cell_Jtx[1])

1;

################################################################

#### DEBUGGING


# function Gridap.Arrays.return_cache(f::PanelRotationField,cellx::AbstractArray{<:VectorValue},panel_id::Int64)
#   A = f.mats[panel_id]
#   x = first(cellx)
#   T = typeof(A⋅x)
#   y = similar(cellx,T)
#   return y # CachedArray(y)
# end

# function Gridap.Arrays.evaluate!(cache,f::PanelRotationField,cellx::AbstractArray{<:VectorValue},panel_id::Int64)
#   # setsize!(cache,size(cellx))
#   y = cache
#   A = f.mats[panel_id]
#   map!(x -> A⋅x, y, cellx)
#   return y
# end

# function Gridap.Arrays.return_cache(f::PanelRotationField,x::VectorValue,panel_id::Int64)
#   y = zero(x)
#   return y
# end

# function Gridap.Arrays.evaluate!(cache,f::PanelRotationField,x::VectorValue,panel_id::Int64)
#   y = cache
#   A = f.mats[panel_id]
#   y = A .⋅ x
#   return y
# end

# """
# if no panel id is provided, assume mats is the correct panel matrix
# """
# function Gridap.Arrays.return_cache(f::PanelRotationField,cellx::AbstractArray{<:VectorValue})
#   A = f.mats
#   x = first(cellx)
#   T = typeof(A⋅x)
#   y = similar(cellx,T)
#   return y # CachedArray(y)
# end

# function Gridap.Arrays.evaluate!(cache,f::PanelRotationField,cellx::AbstractArray{<:VectorValue})
#   # setsize!(cache,size(cellx))
#   y = cache
#   A = f.mats
#   map!(x -> A⋅x, y, cellx)
#   return y
# end







### gradient


# function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:PanelRotationField},cellx::AbstractArray{<:VectorValue},panel_id::Int64)
#   f.object.mats[panel_id]
# end

# function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:PanelRotationField},cellx::AbstractArray{<:VectorValue},panel_id::Int64)
#   w = cache
#   w = f.object.mats[panel_id]
# end

# function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:PanelRotationField},x::VectorValue,panel_id::Int64)
#   f.object.mats[panel_id]
# end

# function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:PanelRotationField},x::VectorValue,panel_id::Int64)
#   w = cache
#   w = f.object.mats[panel_id]
# end

# """
# no panel id -> assume mats is correct
# """
# function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:PanelRotationField},cellx::AbstractArray{<:VectorValue})
#   f.object.mats
# end

# function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:PanelRotationField},cellx::AbstractArray{<:VectorValue})
#   w = cache
#   w = f.object.mats
# end
