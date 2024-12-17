"""
https://github.com/gridap/GridapSolvers.jl/blob/d1cab8c4e829439eb66923a57aebcd4bd69512f5/src/PatchBasedSmoothers/DistributedPatchDecompositions.jl#L83
"""

using GridapGeosciences
using GridapP4est
using PartitionedArrays
using MPI
using Gridap
using GridapDistributed
using DrWatson
using FillArrays
using GridapSolvers
using Gridap.Helpers: @check
import GridapSolvers.MultilevelTools: ModelHierarchyLevel, HierarchicalArray

function adapt_model(ranks,model; refine=true)
  cell_partition=get_cell_gids(model.octree_model.dmodel)
  ref_coarse_flags=map(ranks,partition(cell_partition)) do rank,indices
      flags=zeros(Cint,length(indices))
      if refine
        flags.=refine_flag
      else
        flags.=nothing_flag
      end
  end
  # Gridap.Adaptivity.adapt(model,ref_coarse_flags)
  GridapP4est.adapt(model,ref_coarse_flags)
end

nprocs = (1,1)
ranks = with_mpi() do distribute
  distribute(LinearIndices((prod(nprocs),)))
end

model0 = CubedSphereDiscreteModel(ranks,1;adaptive=true)
model0,_ = adapt_model(ranks,model0)
# model0 = CubedSphereDiscreteModel(4)
# model0 = CartesianDiscreteModel((0,1.0,0.0,1.0), (4,4),isperiodic=(true,true))

face_labeling = get_face_labeling(model0)
topo = Gridap.Geometry.get_grid_topology(model0)
map(local_views(face_labeling),local_views(topo)) do face_labeling, topo
    tag_to_name = face_labeling.tag_to_name
    tag_to_entities = face_labeling.tag_to_entities
    d_to_dface_to_entity = face_labeling.d_to_dface_to_entity
    println(tag_to_name)
    println(tag_to_entities)
    println(d_to_dface_to_entity)
end

    interface_entity = maximum(map(x -> maximum(x;init=0),tag_to_entities)) + 1
    push!(tag_to_entities,[interface_entity])
    push!(tag_to_name,"interface")

    # Interface faces should also be interior
    interior_tag = findfirst(x->(x=="interior"),tag_to_name)
    push!(tag_to_entities[interior_tag],interface_entity)

    # Select interface entities
    boundary_tag = findfirst(x->(x=="boundary"),tag_to_name)
    boundary_entities = tag_to_entities[boundary_tag]


    println(tag_to_name)
    println(tag_to_entities)
# end
Dc = 2
f2c_map = Gridap.Geometry.get_faces(topo,Dc-1,Dc)
    num_cells_around_facet = map(length,f2c_map)
    mx = maximum(num_cells_around_facet)
    for (f,nf) in enumerate(num_cells_around_facet)
      is_boundary = (d_to_dface_to_entity[Dc][f] ∈ boundary_entities)
      if !is_boundary && (nf != mx)
        d_to_dface_to_entity[Dc][f] = interface_entity
      end
    end



model0 = CubedSphereDiscreteModel(4)
PatchDecomposition(model0)
