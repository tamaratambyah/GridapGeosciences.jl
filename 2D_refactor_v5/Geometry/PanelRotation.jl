"""
PanelRotation
- PanelRotationMap
- PanelRotationField

Applies the roation provided by mats to x
- `mats` can be either a single TensorValue, or an array of TensorValues.
- if `panel_id` is not provided, assume `mats` is the correct panel matrix

The same functionality is provided for a GridapMap and GridapField
- the GridapField is used to map cell coordinates and cell maps
    * the gradient is implemented
- the GridapMap is generally not used, but is provided regardless for testing
purposes

This interface can also be used as the inverse if mats^{-1} is the provided input

"""

"""
PanelRotationField
"""

struct PanelRotationField{A} <: Gridap.Fields.Field
  mats::A
end


function Gridap.Arrays.return_cache(f::PanelRotationField,cellx::AbstractArray{<:VectorValue},panel_id::Int64)
  A = f.mats[panel_id]
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationField,cellx::AbstractArray{<:VectorValue},panel_id::Int64)
  y = cache
  A = f.mats[panel_id]
  map!(x -> A⋅x, y, cellx)
  return y
end

function Gridap.Arrays.return_cache(f::PanelRotationField,x::VectorValue,panel_id::Int64)
  y = zero(x)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationField,x::VectorValue,panel_id::Int64)
  y = cache
  A = f.mats[panel_id]
  y = A .⋅ x
  return y
end

"""
if no panel id is provided, assume mats is the correct panel matrix
"""
function Gridap.Arrays.return_cache(f::PanelRotationField,cellx::AbstractArray{<:VectorValue})
  A = f.mats
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationField,cellx::AbstractArray{<:VectorValue})
  y = cache
  A = f.mats
  map!(x -> A⋅x, y, cellx)
  return y
end

"""
if no panel id is provided, assume mats is the correct panel matrix
"""
function Gridap.Arrays.return_cache(f::PanelRotationField,x::VectorValue)
  y = zero(x)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationField,x::VectorValue)
  y = cache
  A = f.mats
  y = A .⋅ x
  return y
end


## gradient: is just the coefficient matrix (mats)
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)

function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:PanelRotationField},cellx::AbstractArray{<:VectorValue},panel_id::Int64)
  A = f.object.mats[panel_id]
  T = typeof(transpose(A))
  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:PanelRotationField},cellx::AbstractArray{<:VectorValue},panel_id::Int64)

  y = cache
  A = f.object.mats[panel_id]
  map!(x -> transpose(A), y, cellx)
  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:PanelRotationField},x::VectorValue,panel_id::Int64)
  transpose( f.object.mats[panel_id] )
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:PanelRotationField},x::VectorValue,panel_id::Int64)
  w = cache
  w = f.object.mats[panel_id]
  return transpose(w)
end


"""
if no panel id is provided, assume mats is correct
"""

function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:PanelRotationField},cellx::AbstractArray{<:VectorValue})
  println("getting cache")
  A = f.object.mats
  T = typeof(transpose(A))
  y = similar(cellx,T,size(cellx))
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:PanelRotationField},cellx::AbstractArray{<:VectorValue})
  c, = cache
  setsize!(c,size(cellx))
  y = c.array
  A = f.object.mats
  # map!(x -> transpose(A), y, cellx)
  fill!(y,transpose(A))
  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:PanelRotationField},x::VectorValue)
  f.object.mats
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:PanelRotationField},x::VectorValue)
  w = cache
  w = f.object.mats
  return transpose(w)
end


"""
PanelRotationMap
"""
struct PanelRotationMap{A} <: Map
  mats::A
end

Rp1PanelMap() = PanelRotationMap(rp1_3D)
R1pPanelMap() = PanelRotationMap(r1p_3D)

function Gridap.Arrays.return_cache(f::PanelRotationMap,cellx::AbstractArray{<:VectorValue},panel_id::Int64)
  A = f.mats[panel_id]
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,cellx::AbstractArray{<:VectorValue},panel_id::Int64)
  y = cache
  A = f.mats[panel_id]
  map!(x -> A⋅x, y, cellx)
  return y
end

function Gridap.Arrays.return_cache(f::PanelRotationMap,x::VectorValue,panel_id::Int64)
  y = zero(x)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,x::VectorValue,panel_id::Int64)
  y = cache
  A = f.mats[panel_id]
  y = A .⋅ x
  return y
end

"""
if no panel id is provided, assume mats is the correct panel matrix
"""
function Gridap.Arrays.return_cache(f::PanelRotationMap,cellx::AbstractArray{<:VectorValue})
  A = f.mats
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,cellx::AbstractArray{<:VectorValue})
  y = cache
  A = f.mats
  map!(x -> A⋅x, y, cellx)
  return y
end

"""
if no panel id is provided, assume mats is the correct panel matrix
"""
function Gridap.Arrays.return_cache(f::PanelRotationMap,x::VectorValue)
  y = zero(x)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,x::VectorValue)
  y = cache
  A = f.mats
  y = A .⋅ x
  return y
end



"""
panel_rotations

returns the TensorValue matrices required to map points from panel 1 <-> panel i
"""


function panel_rotations(d::Int)
  if d == 2
    rotate_panel_p_to_1, rotate_panel_1_to_p = _panel_rotation_matrices_2D()
  elseif d == 3
    rotate_panel_p_to_1, rotate_panel_1_to_p = _panel_rotation_matrices_3D()
  end
  rp1 = map(TensorValue,rotate_panel_p_to_1)
  r1p = map(TensorValue,rotate_panel_1_to_p)
  return rp1, r1p
end

function _panel_rotation_matrices_3D()
  # rotation about X axis by 3π/2
  Rx = [1.0 0.0 0.0
        0.0 0.0 1.0
        0.0 -1.0 0.0]

  # rotation about Y axis by 3π/2
  Ry = [0.0 0.0 -1.0
        0.0 1.0 0.0
        1.0 0.0 0.0]

  # rotation about Z axis by π/2
  Rz = [0.0 -1.0 0.0
        1.0 0.0 0.0
        0.0 0.0 1.0]


  """ Panel 1: front  """
  A_11 = Matrix{Float64}(I,3,3)


  """ Panel 2: top"""
  A_12 = Ry # panel 1 -> 2: rotation about Y axis by 3π/2
  A_21 = inv(A_12) # panel 2 -> 1

  """ Panel 3: right """
  A_23 = Rx # panel 2 -> 3: rotation about X axis by 3π/2
  A_13 = A_23*A_12 # panel 1-> 3: panel 2 -> 3 ⋅ panel 1 -> 2
  A_31 = inv(A_13) # panel 3 -> 1

  """ Panel 4: back """
  A_34 = Rz # panel 3 -> 4: rotation about Z axis by π/2
  A_14 = A_34*A_13 # panel 1 -> 4: panel 3 -> 4 ⋅ panel 1 -> 3
  A_41 = inv(A_14) # panel 4 -> 1

  """ Panel 5: bottom"""
  A_45 = Ry # panel 4 -> 5: rotation about Y axis by 3π/2
  A_15 = A_45*A_14 # panel 1 -> 5: panel 4 -> 5 ⋅ panel 1 -> 4
  A_51 = inv(A_15)

  """ Panel 6: left"""
  A_56 = Rx # panel 5 -> 6: rotation about X axis by 3π/2
  A_16 = A_56*A_15 # panel 1 -> 6: panel 5 -> 6 ⋅ panel 1 -> 5
  A_61 = inv(A_16)

  rotate_panel_p_to_1 = [A_11, A_21, A_31, A_41, A_51, A_61]
  rotate_panel_1_to_p = [A_11, A_12, A_13, A_14, A_15, A_16]

  return rotate_panel_p_to_1, rotate_panel_1_to_p
end

function _panel_rotation_matrices_2D()
  # rotation by π/2
  R90 = [0.0 -1.0
        1.0 0.0]

  # rotation by π
  R180 = [-1.0 0.0
        0.0 -1.0]


  """ Panel 1: front  """
  A_11 = Matrix{Float64}(I,2,2)  # panel 1 -> 1: no rotation


  """ Panel 2: top"""
  A_12 = Matrix{Float64}(I,2,2) # panel 1 -> 2: no rotation
  A_21 = inv(A_12) # panel 2 -> 1

  """ Panel 3: right """
  A_13 = R90 # panel 1-> 3: rotate 90
  A_31 = inv(A_13) # panel 3 -> 1

  """ Panel 4: back """
  A_14 = R180 # panel 1 -> 4: rotate 180
  A_41 = inv(A_14) # panel 4 -> 1

  """ Panel 5: bottom"""
  A_15 = R180 # panel 1 -> 5: rotate 180
  A_51 = inv(A_15)

  """ Panel 6: left"""
  A_16 = R90 # panel 1 -> 6: rotate 90
  A_61 = inv(A_16)

  rotate_panel_p_to_1 = [A_11, A_21, A_31, A_41, A_51, A_61]
  rotate_panel_1_to_p = [A_11, A_12, A_13, A_14, A_15, A_16]

  return rotate_panel_p_to_1, rotate_panel_1_to_p
end

const rp1_3D, r1p_3D = panel_rotations(3)
const rp1_2D, r1p_2D = panel_rotations(2)
