function bump_matrics()
  _B = [0.0 0.0
        1.0 0.0
        0.0 1.0]
  B = TensorValue(_B)
  A = transpose(B)
  b = Point(1.0,0.0,0.0)
  return A, B, b
end

# function bump_panel1_2D_to_3D(x::Point{2})
#   A,B,b = bump_matrics()
#   B⋅x + b
# end

# function bump_panel1_3D_to_2D(x::Point{3})
#   A,B,b = bump_matrics()
#   A⋅x
# end


# A,B,b = bump_matrics()

# bump_panel1_2D_to_3D(Point(1,1))
# bump_panel1_3D_to_2D(Point(1,1,1))
