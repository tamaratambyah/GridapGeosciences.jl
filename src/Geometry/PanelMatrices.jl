function cube_to_αβ(p::Int)
  if p == 1
    A = [0 1 0
        0 0 1]
  elseif p == 2
    A = [0 1 0
        -1 0 0]
  elseif p == 3
    A = [-1 0 0
          0 0 1]
  elseif p == 4
    A = [0 0 1
        0 1 0]
  elseif p == 5
    A = [-1 0 0
        0 1 0]
  elseif p == 6
    A = [0 0 1
        -1 0 0]
  end

  return TensorValue(A)

end

function cube_face_values(p::Int;a::Float64=CUBE_HALF_EDGE)
  if p == 1 # X = 1 panel (front)
    b = VectorValue(1.0, 0.0, 0.0)
  elseif p == 2 # Z = 1 panel (top)
    b = VectorValue(0.0, 0.0, 1.0)
  elseif p == 3 # Y = 1 panel (right)
    b = VectorValue(0.0, 1.0, 0.0)
  elseif p == 4 # X = -1 panel (back)
    b = VectorValue(-1.0, 0.0, 0.0)
  elseif p == 5 # Z = -1 panel (bottom)
    b = VectorValue(0.0, 0.0, -1.0)
  elseif p == 6 # Y = -1 panel (left)
    b = VectorValue(0.0, -1.0, 0.0)
  end

  return a*b

end

function rotation_mats(p::Int)
  if p == 1
    A = [1 0 0
         0 1 0
         0 0 1]
  elseif p == 2
    A = [0 0 -1
         0 1 0
         1 0 0]
  elseif p == 3
    A = [0 -1 0
         1 0 0
         0 0 1]
  elseif p == 4
    A = [-1 0 0
          0 0 1
          0 1 0]
  elseif p == 5
    A = [0 -1 0
         0 0 1
        -1 0 0]
  elseif p == 6
    A = [0 0 -1
         -1 0 0
         0 1 0]
  end
  TensorValue(A)
end

R1p = map(p->rotation_mats(p),collect(1:6))

A_cube2panel = map(p->cube_to_αβ(p),collect(1:6))
A_panel2cube = map(A->transpose(A),A_cube2panel)
b_panel2cube = map(p->cube_face_values(p),collect(1:6))
