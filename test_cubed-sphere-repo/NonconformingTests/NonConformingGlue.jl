using GridapP4est
using Gridap
using PartitionedArrays
using GridapDistributed
using MPI
using Gridap.FESpaces
using FillArrays
using Test
using DrWatson
using Gridap.Geometry

MPI.Init()
np = MPI.Comm_size(MPI.COMM_WORLD)
ranks = distribute_with_mpi(LinearIndices((np,)))

include("refinement_helpers.jl")

coarse_model = CartesianDiscreteModel((0,1,0,1),(3,3),isperiodic=(true,true))
dmodel = OctreeDistributedDiscreteModel(ranks,coarse_model)#,2)

ref_coarse_flags = middle_refinement(dmodel)
# ref_coarse_flags = boundary_refinement(dmodel)
fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);



lmodel = fmodel.dmodel.models.item
topo = get_grid_topology(lmodel)
Dc = num_cell_dims(topo)
Gridap.Geometry.get_faces(topo,Dc,0) # nodes
Gridap.Geometry.get_faces(topo,Dc,1) # faces

non_conforming_glue = fmodel.non_conforming_glue
gridap_cell_faces = map(local_views(fmodel.dmodel)) do model
        topo = Gridap.Geometry.get_grid_topology(model)
        Tuple(Gridap.Geometry.get_faces(topo, Dc, d) for d = 0:Dc-1)
    end


num_regular_faces = map(non_conforming_glue) do ncglue
    println("num_regular_faces=$(Tuple(ncglue.num_regular_faces[d] for d = 1:Dc))")
    Tuple(ncglue.num_regular_faces[d] for d = 1:Dc)
end

num_hanging_faces = map(non_conforming_glue) do ncglue
    println("num_hanging_faces=$(Tuple(ncglue.num_hanging_faces[d] for d = 1:Dc))")
    Tuple(ncglue.num_hanging_faces[d] for d = 1:Dc)
end

hanging_faces_glue = map(non_conforming_glue) do ncglue
    Tuple(ncglue.hanging_faces_glue[d] for d = 1:Dc)
end

hanging_faces_to_cell = map(non_conforming_glue) do ncglue
    Tuple(ncglue.hanging_faces_to_cell[d] for d = 1:Dc)
end

hanging_faces_to_lface = map(non_conforming_glue) do ncglue
    Tuple(ncglue.hanging_faces_to_lface[d] for d = 1:Dc)
end

owner_faces_pindex=map(non_conforming_glue) do ncglue
    Tuple(ncglue.owner_faces_pindex[d] for d = 1:Dc-1)
end

owner_faces_lids=map(non_conforming_glue) do ncglue
    Tuple(ncglue.owner_faces_lids[d] for d = 1:Dc-1)
end
