"""
pinv

left pseudo inverse of general D1xD2 TensorValue, J,  defined as
J^† = (J^T J)^{-1} J^T such that J^† J = I
"""
function pinv(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  @check D1 !== D2
  Jt = transpose(J)
  inv(Jt⋅J)⋅Jt
end
