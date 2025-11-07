using Gridap
using GridapDistributed
using GridapGeosciences
using MPI
using PartitionedArrays
using DrWatson
include("../convergence_tools.jl")

MPI.Init()

nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
# nprocs = 2
# ranks = with_debug() do distribute
#   distribute(LinearIndices((nprocs,)))
# end

models =  get_distributed_refined_models(ranks,nprocs,1,true)

panel_model = models[1]
i_am_main(ranks) && println(num_cells(panel_model))
Ω_panel = Triangulation(panel_model)
panel_ids = get_panel_ids(panel_model)

meas_cf_skel = CellField(_sqrtg,Ω_panel)
Λ = SkeletonTriangulation(with_ghost,panel_model)

degree = 0
dΛ = Measure(Λ,degree)
dΩ = Measure(Ω_panel,degree)

function f(rank, meas_cf_skel_loc, dΛ_loc)
  println(typeof(meas_cf_skel_loc))
  quad = dΛ_loc.quad
  field = Gridap.CellData.change_domain(meas_cf_skel_loc,quad.trian,quad.data_domain_style)
  if rank==1
    cell_array_fields = Gridap.CellData.get_data(field)
    i_am_main(rank) && println(typeof(cell_array_fields[end]))
    i_am_main(rank) && println(collect(field(get_cell_points(quad))))
  end
end

map(f, ranks, local_views(meas_cf_skel), local_views(dΛ))
