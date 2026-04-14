
abstract type TSWODEOperator{T} <: ODEOperator{T} end
ODEOperatorType(::TSWODEOperator{T}) where {T} = T
ODEOperatorType(::Type{<:TSWODEOperator{T}}) where {T} = T
