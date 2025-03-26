"""
bump_matrices

returns the matrices required to bump points 2D <-> 3D on panel 1.
note, it is assumed the cube faces are [-a,a]^2 so that panel 1 has X = a
"""
function bump_matrics(a::Float64)
  A = [0.0 1.0 0.0
       0.0 0.0 1.0]
  B = [0.0 0.0
       1.0 0.0
       0.0 1.0]
  b = a .* [1.0
                0.0
                0.0]
  return A, B, b
end

"""
panel_rotations

returns the matrices required to map points from panel 1 <-> panel i
"""
function panel_rotations()
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
