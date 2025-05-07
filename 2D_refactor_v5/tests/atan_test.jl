x = 1.0
y = 1.0

gle = my_atan(-0.0,-0.0)

function my_atan(x,y)


  if x == 0.0 && y == 0.0
    println("centre")
    return angle = 0.0
  end


  # x axis
  if y == 0.0
    println("x axis")
    if x > 0.0
      return angle = 0.0
    elseif x < 0.0
      return angle = π
    end
  end

  # y axis
  if x == 0.0
    println("y axis")
    if y > 0.0
      return angle = π/2
    elseif y < 0.0
      return angle = 3*π/2
    end
  end



  ref_angle = atan(abs(y)/abs(x))
  if x > 0.0 && y > 0.0 # quad 1
    println("quad 1")
    angle = ref_angle
  elseif x < 0.0 && y > 0.0 # quad 2
    println("quad 2")
    angle = π -  ref_angle
  elseif x < 0.0 && y < 0.0 # quad 3
    println("quad 3")
    angle = π +  ref_angle
  elseif x > 0.0 && y < 0.0 # quad 4
    println("quad 4")
    angle = 2*π -  ref_angle
  end

  return angle
end
