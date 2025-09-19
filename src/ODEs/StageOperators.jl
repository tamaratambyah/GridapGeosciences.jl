
#######################
# LinearStageOperator #
#######################


function LinearStageOperator(
  odeop::DAEODEOperator, odeopcache,
  tx::Real, usx::Tuple{Vararg{AbstractVector}},
  ws::Tuple{Vararg{Real}},
  J::AbstractMatrix, r::AbstractVector, reuse::Bool, sysslvrcache
)
  Gridap.ODEs.residual!(r, odeop, tx, usx, odeopcache)

  if isnothing(sysslvrcache) || !reuse
    Gridap.ODEs.jacobian!(J, odeop, tx, usx, ws, odeopcache)
  end

  LinearStageOperator(J, r, reuse)
end
