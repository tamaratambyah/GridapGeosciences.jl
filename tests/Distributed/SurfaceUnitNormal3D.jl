using DrWatson
using Gridap
using GridapDistributed
using GridapP4est
using GridapGeosciences
using GridapGeosciences.Distributed

using MPI
using PartitionedArrays
using Test


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

function test_normal_unit_vector(panel_model::GridapDistributed.DistributedDiscreteModel{3},dir,return_vtk=false)
  ranks = get_ranks(panel_model)

  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v")

  Ω_panel =  Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)
  dΩ = Measure(Ω_panel,4)

  ## the normal in parametric space (γ,α,β) is (1,0,0)
  n3D_panel = CellField(VectorValue(1.0,0.0,0.0),Ω_panel)
  J_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)
  inv_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)

  ## map the normal from parametric space -> ambient space
  _n_mapped = J_cf ⋅ (inv_cf  ⋅ n3D_panel )
  ff = Operation(sqrt)(  n3D_panel   ⋅ (inv_cf⋅ n3D_panel )  )
  n_mapped = _n_mapped/ff

  ## the unit surface normal is given by the position vector
  vX = panel_to_cartesian(normal_vec)
  norm_vec_cf = panelwise_cellfield(vX,Ω_panel,panel_ids)

  e = l2(norm_vec_cf-n_mapped,dΩ)
  i_am_main(ranks) && println(e)
  @test e < 1e-12

  if return_vtk
    lvl = nref(panel_model)
    cell_geo_map = geo_map_func(Ω_panel)
    panel_cfs = [ n_mapped,norm_vec_cf,norm_vec_cf-n_mapped]
    labels = ["n_mapped", "n_vec", "diff"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

end

include("../convergence_tools.jl")
dir = datadir("NormalVector3D")
!isdir(dir) && mkdir(dir)

models = get_3D_octree_refined_models(ranks,2)
return_vtk = false

for panel_model in models
  test_normal_unit_vector(panel_model,dir,return_vtk)
end

# ./mpiexecjl.sh -np 2 julia --project=. $PWD/tests/SurfaceUnitNormal3D.jl
