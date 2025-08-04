using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields, Gridap.TensorValues

################################################################################
#### I think the following extra overloads should eventually go to Gridap's core
#### They are needed to have the product rule working whenever there is the
#### specific combination of Field and ArrayBlock{<:Field} that is triggered below
################################################################################
function Gridap.Fields.return_value(k::Broadcasting{<:Operation},
                                    f1::Field,f2::ArrayBlock{A,N},g1::Field,g2::ArrayBlock{B,N}) where {A,B,N}
  @assert length(f2.array) == length(g2.array)
  @assert f2.touched == g2.touched

  f2i = Gridap.Fields.testitem(f2)
  g2i = Gridap.Fields.testitem(g2)
  f1f2ig1g2i = Gridap.Fields.return_value(k,f1,f2i,g1,g2i)
  o = Array{typeof(f1f2ig1g2i),N}(undef,size(f2.array))
  for i in eachindex(f2.array)
    if f2.touched[i]
      o[i] = Gridap.Fields.return_value(k,f1,f2.array[i],g1,g2.array[i])
    end
  end
  ArrayBlock(o,f2.touched)
end

function Gridap.Fields.return_cache(k::Broadcasting{<:Operation},
                                    f1::Field,f2::ArrayBlock{A,N},g1::Field,g2::ArrayBlock{B,N}) where {A,B,N}
  @assert length(f2.array) == length(g2.array)
  @assert f2.touched == g2.touched

  f2i = Gridap.Fields.testitem(f2)
  g2i = Gridap.Fields.testitem(g2)

  cf1f2ig1g2i = Gridap.Fields.return_cache(k,f1,f2i,g1,g2i)
  f1f2ig1g2i = Gridap.Fields.evaluate!(cf1f2ig1g2i,k, f1,f2i,g1,g2i)

  l = Array{typeof(cf1f2ig1g2i),N}(undef,size(f2.array))
  o = Array{typeof(f1f2ig1g2i),N}(undef,size(f2.array))
  for i in eachindex(f2.array)
    if f2.touched[i]
      l[i] = return_cache(k,f1,f2.array[i],g2,g2.array[i])
    end
  end
  ArrayBlock(o,f2.touched),l
end

function Gridap.Fields.evaluate!(cache,k::Broadcasting{<:Operation},
       f1::Field,f2::ArrayBlock{A,N},g1::Field,g2::ArrayBlock{B,N}) where {A,B,N}

  o,l = cache
  @assert o.touched == f2.touched
  for i in eachindex(f2.array)
     if f2.touched[i]
       o.array[i] = evaluate!(l[i],k,f1, f2.array[i],g1,g2.array[i])
     end
  end
  o
end
################################################################################
#### End of extra overloads
################################################################################
