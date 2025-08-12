function cube_to_αβ(p)
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

function rotation_mats(p)
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
