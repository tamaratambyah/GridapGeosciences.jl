function Jfactor(x)
  α,β = x
  RADIUS/( (1 + (tan(α))^2 + (tan(β))^2 )^(3/2) )
  1.0
end

function dXdalpha(x)
  α,β = x
  Jfactor(x)* ( -1.0*tan(α)*(sec(α))^2 )
end

function dXdbeta(x)
  α,β = x
  Jfactor(x)* ( -1.0*tan(β)*(sec(β))^2 )
end

function dYdalpha(x)
  α,β = x
  Jfactor(x)* ( (sec(α))^2*(sec(β))^2 )
end

function dYdbeta(x)
  α,β = x
  Jfactor(x)* ( -1.0*tan(α)*tan(β)*(sec(β))^2 )
end

function dZdalpha(x)
  α,β = x
  Jfactor(x)* ( -1.0*tan(α)*tan(β)*(sec(α))^2 )
end

function dZdbeta(x)
  α,β = x
  Jfactor(x)* ( (sec(α))^2*(sec(β))^2 )
end

function Jpanel1(x)
  TensorValue{3,2}( dXdalpha(x), dYdalpha(x), dZdalpha(x),
                    dXdbeta(x), dYdbeta(x), dZdbeta(x)
  )
end

function JTpanel1(x)
  TensorValue{2,3}( dXdalpha(x), dXdbeta(x),
                    dYdalpha(x), dYdbeta(x),
                    dZdalpha(x), dZdbeta(x)
  )
end


function Jpanel2(x)
  TensorValue{3,2}( -1.0*dZdalpha(x), dYdalpha(x), dXdalpha(x),
                    -1.0*dZdbeta(x), dYdbeta(x), dXdbeta(x)
  )
end

function Jpanel3(x)
  TensorValue{3,2}( -1.0*dYdalpha(x), dXdalpha(x), dZdalpha(x),
                    -1.0*dYdbeta(x), dXdbeta(x), dZdbeta(x)
  )
end

function Jpanel4(x)
  TensorValue{3,2}( -1.0*dXdalpha(x), dZdalpha(x), dYdalpha(x),
                    -1.0*dXdbeta(x), dZdbeta(x), dYdbeta(x)
  )
end

function Jpanel5(x)
  TensorValue{3,2}( -1.0*dYdalpha(x), dZdalpha(x), -1.0*dXdalpha(x),
                    -1.0*dYdbeta(x), dZdbeta(x), -1.0*dXdbeta(x)
  )
end

function Jpanel6(x)
  TensorValue{3,2}( -1.0*dZdalpha(x), -1.0*dXdalpha(x), dYdalpha(x),
                    -1.0*dZdbeta(x), -1.0*dXdbeta(x), dYdbeta(x)
  )
end

##################################################################################
### Panel 2
function Jfactor2(x)
  α,β = x
  RADIUS/( (1 + (tan(α))^2 + (cot(β))^2 )^(3/2) )
end

function dXdalpha2(x)
  α,β = x
  Jfactor2(x)* ( -1.0*tan(α)*(sec(α))^2 )
end

function dXdbeta2(x)
  α,β = x
  Jfactor2(x)* ( cot(β)*(csc(β))^2 )
end

function dYdalpha2(x)
  α,β = x
  Jfactor2(x)* ( (sec(α))^2*(csc(β))^2 )
end

function dYdbeta2(x)
  α,β = x
  Jfactor2(x)* ( tan(α)*cot(β)*(csc(β))^2 )
end

function dZdalpha2(x)
  α,β = x
  Jfactor2(x)* ( tan(α)*cot(β)*(sec(α))^2 )
end

function dZdbeta2(x)
  α,β = x
  Jfactor2(x)* ( (sec(α))^2*(csc(β))^2 )
end


function Jpanel2_new(x)
  TensorValue{3,2}( -1.0*dZdalpha2(x), dYdalpha2(x), dXdalpha2(x),
                    -1.0*dZdbeta2(x), dYdbeta2(x), dXdbeta2(x)
  )
end
