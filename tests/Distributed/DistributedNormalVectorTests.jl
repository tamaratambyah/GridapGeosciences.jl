"""
Test the pushforward of normal vectors, pullback of area form, and continuity
of metric on skeleton mesh
"""

module DistributedNormalTests

using Gridap
using GridapGeosciences
using GridapDistributed
using Test

include("../convergence_tools.jl")

################################################################################
#### Test unit normal vectors
################################################################################
myisless(b::Gridap.TensorValues.MultiValue,a::Number) = all(Gridap.TensorValues.isless.(b.data,a))

function test_debug_vector_equality(out,tol=1e-14)
  map(out) do o
    @test all( lazy_map(x-> all(myisless.(x,tol)), o))
  end
end

function test_debug_equality(out,tol=1e-14)
  map(out) do o
    @test all( lazy_map(x-> all(isless.(x,tol)), o))
  end
end

function main(distribute,nprocs)

  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 2
  dmodels = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)
  panel_model = dmodels[2]

  Ω_panel = Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)
  Λ = SkeletonTriangulation(with_ghost,panel_model)
  n_Λ = get_normal_vector(Λ)
  pts = get_cell_points(Λ)
  ##############################################################################
  ## Face normals: skeleton trian
  ##############################################################################

  # Method 1: Use gridap machinary
  n = pushforward_normal(Λ)
  out = (n.plus+n.minus)(pts)
  test_debug_vector_equality(out)

  # Method 2: Santi's formula
  cell_geo_map = geo_map_func(panel_ids)
  n = pushforward_normal(Λ,cell_geo_map)
  out = (n.plus+n.minus)(pts)
  test_debug_vector_equality(out)

  ##############################################################################
  ## DG tests
  ### check sqrt(g) is continuous across skeleton
  ### check |Jg^-1 n| - pullback of area form
  ##############################################################################
  meas_cf = panelwise_cellfield(_sqrtg,Λ)
  out = (meas_cf.plus-meas_cf.minus)(pts)
  test_debug_equality(out)

  area_form_cf = pullback_area_form(Λ)
  out = (area_form_cf.plus-area_form_cf.minus)(pts)
  test_debug_equality(out)

  ##############################################################################
  ## Advection tests
  ### check abs(v⋅n.plus) = abs(v⋅n.minus)
  ##############################################################################
  vecX(XYZ) = VectorValue(-XYZ[2],XYZ[3],0.0)
  vX = panel_to_cartesian(tangent_vec(vecX))

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)

  _vel = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  vel = interpolate(_vel,U)

  diff_cf = (abs((vel⋅ n_Λ).minus) .- abs((vel⋅ n_Λ).plus))(pts)
  test_debug_equality(diff_cf)

end


end # module
