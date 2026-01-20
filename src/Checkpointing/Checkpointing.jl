module Checkpointing

using PartitionedArrays
using GridapDistributed
using JLD2
using Gridap.FESpaces
using GridapDistributed

function psave(dir::AbstractString, xh::FEFunction)
  psave(dir,xh.free_values)
end

function psave(dir::AbstractString, x::AbstractArray)
  !isdir(dir) && mkdir(dir)
  arr = to_local_storage(x)
  filename = joinpath(dir,basename(dir)*".jld2")
  save_object(filename,arr)
end

function pload(dir::AbstractString,ranks::AbstractArray{Bool})
  filename = joinpath(dir,basename(dir)*".jld2")
  arr = load_object(filename)
  return from_local_storage(arr)
end



"""
    psave(filename::AbstractString, x)

Save a partitioned object `x` to a directory `dir` as
a set of JLD2 files corresponding to each part.
"""
function psave(dir::AbstractString,
  xh::Union{GridapDistributed.DistributedCellField,GridapDistributed.DistributedMultiFieldCellField})
  psave(dir,xh.metadata.free_values)
end

function psave(dir::AbstractString, x::Union{PVector,PSparseMatrix})
  ranks = get_parts(x)
  i_am_main(ranks) && mkpath(dir)
  arr = to_local_storage(x)
  i_am_main(ranks) && println("saving solution")

  # ensure no MPI task tries to generate the file before the main MPI task has
  # created the folder
  PartitionedArrays.barrier(ranks)

  # write each file
  map(ranks,arr) do id, arr
    println(id)
    filename = joinpath(dir,basename(dir)*"_$id.jld2")
    # save_object(filename,arr)
    # jldsave(filename; arr)
    jldopen(filename, "w+") do f
      f["$id"] = arr
    end
  end
  # PartitionedArrays.barrier(ranks)
end

"""
    pload(dir::AbstractString, ranks::AbstractArray{<:Integer})

Load a partitioned object stored in a set of JLD2 files in directory `dir`
indexed by MPI ranks `ranks`.
"""
function pload(dir::AbstractString, ranks::AbstractArray{<:Integer})
  arr = map(ranks) do id
    filename = joinpath(dir,basename(dir)*"_$id.jld2")
    # load_object(filename)
    _A = jldopen(filename, "r") do _f
      _f["$id"]
    end
    return _A
  end
  y = from_local_storage(arr)
  GridapDistributed.DistributedFEFunctionData(y).free_values
end

"""
    pload!(dir::AbstractString, x)

Load a partitioned object stored in a set of JLD2 files in directory `dir`
and copy contents to the equivilent object `x`.
"""
function pload!(dir::AbstractString, x)
  ranks = get_parts(x)
  y = pload(dir,ranks)
  copyto!(x,y)
  return x
end

## Handle storage of values and indices for `PVector` and PSparseMatrix

function to_local_storage(x)
  x
end

function from_local_storage(x)
  x
end

struct PVectorLocalStorage{A,B}
  values::A
  indices::B
end

function from_local_storage(x::AbstractArray{<:PVectorLocalStorage})
  values, indices = map(x) do x
    x.values, x.indices
  end |> tuple_of_arrays
  return PVector(values,indices)
end

function to_local_storage(x::PVector)
  map(PVectorLocalStorage,partition(x),partition(axes(x,1)))
end

struct PSparseMatrixLocalStorage{A,B,C}
  values::A
  rows::B
  cols::C
end

function from_local_storage(x::AbstractArray{<:PSparseMatrixLocalStorage})
  values, rows, cols = map(x) do x
    x.values, x.rows, x.cols
  end |> tuple_of_arrays
  return PSparseMatrix(values,rows,cols)
end

function to_local_storage(x::PSparseMatrix)
  map(PSparseMatrixLocalStorage,partition(x),partition(axes(x,1)),partition(axes(x,2)))
end


end
