module NormalTests

using DrWatson
using Gridap
using GridapGeosciences
using Test
using Gridap.Helpers

include("../convergence_tools.jl")

################################################################################
#### Test unit normal vectors
################################################################################

function test_normal_unit_vector(panel_model,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  Ω_panel = Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)
  dΩ = Measure(Ω_panel,4)

  vX = panel_to_cartesian(normal_vec)
  norm_vec_cf = panelwise_cellfield(vX,Ω_panel,panel_ids)
  norm_vec_from_basis_cf = panelwise_cellfield(normal_vector_from_basis,Ω_panel,panel_ids)

  @check l2(norm_vec_cf-norm_vec_from_basis_cf,dΩ) < 1e-12

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = geo_map_func(panel_ids)
    panel_cfs = [ norm_vec_cf,norm_vec_from_basis_cf]
    labels = ["normal", "n"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end
end


n_ref_lvls = 4
models  = get_refined_models(n_ref_lvls)
return_vtk = false
dir = datadir("NormalTests")
(!isdir(dir) && return_vtk) && mkdir(dir)

################################################################################
## Unit normal to surface: k = a₁ × a₂
################################################################################
for panel_model in models
  test_normal_unit_vector(panel_model)
end

################################################################################
## Face normals -- for 1 model
################################################################################
panel_model = models[4]
panel_ids = get_panel_ids(panel_model)
cell_geo_map = geo_map_func(panel_ids)

################################################################################
## Face normals: boundary trian
################################################################################
trian = BoundaryTriangulation(panel_model,1)
pts = get_cell_points(trian)

# get face normals in 3D
n_3D = pushforward_normal(trian,cell_geo_map)

# push forward 2D chart normals to ambient space
n_mapped = pushforward_normal(trian)

# test equality
@test all(n_mapped(pts) .≈ n_3D(pts))

################################################################################
## Face normals: skeleton trian
################################################################################
trian = SkeletonTriangulation(panel_model)
pts = get_cell_points(trian)

# regular 2D normal in chart
n_2D = get_normal_vector(trian)

# # get face normals in 3D
# cell_vectors = get_facet_normal(trian,cell_geo_map)
# n_3D = get_normal_vector(trian,cell_vectors)
n_3D = pushforward_normal(trian,cell_geo_map)

# push forward 2D chart normals to ambient space
n_mapped = pushforward_normal(trian)

# test equality of plus and minus side
@test all(n_mapped.plus(pts) .≈ n_3D.plus(pts))
@test all(n_mapped.minus(pts) .≈ n_3D.minus(pts))

## plot normals on skeleton
if return_vtk
  panel_cfs = [n_3D.plus, n_3D.minus, n_3D.minus+n_3D.plus,
              n_2D.plus, n_2D.minus, n_2D.minus+n_2D.plus]
  labels = ["amb_n_plus", "amb_n_minus", "amb_n_total",
            "chart_n_plus", "chart_n_minus", "chart_n_total"]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)

  skel_panel_ids = get_panel_ids(trian)
  skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
  writevtk(trian,dir*"/ambient_model_skeleton",cellfields=cellfields,append=false,geo_map=skel_geo_map)
end
################################################################################
## DG tests
### check sqrt(g) is continuous across skeleton
### show g, g^-1, J not continuous
### check |Jg^-1 n| - pullback of area form
################################################################################
Λ = SkeletonTriangulation(panel_model)
pts = get_cell_points(Λ)
n_Λ = get_normal_vector(Λ)
meas_cf = panelwise_cellfield(sqrtg,Λ)
inv_metric_cf = panelwise_cellfield(inv_metric,Λ)
jac_cf = panelwise_cellfield(forward_jacobian,Λ)
area_form_cf = pullback_area_form(Λ)

# test equality of plus and minus side of sqrt(g)
@test all(meas_cf.minus(pts) .≈ meas_cf.plus(pts))

# test equality of plus and minus side of |Jg^-1 n|
@test all(area_form_cf.plus(pts) .≈ area_form_cf.minus(pts))

# test inequality of plus and minus side for g^-1 and J
@test !all(inv_metric_cf.plus(pts) .≈ inv_metric_cf.minus(pts))
@test !all(jac_cf.minus(pts) .≈ jac_cf.plus(pts))

if return_vtk
  panel_cfs = [meas_cf.plus, meas_cf.minus, meas_cf.minus-meas_cf.plus,
              jac_cf.plus, jac_cf.minus, jac_cf.minus-jac_cf.plus,
              area_form_cf.plus, area_form_cf.minus, area_form_cf.plus-area_form_cf.minus]
  labels = ["g_plus", "g_minus", "g_diff", "jac_plus", "jac_minus", "jac_diff",
            "a_plus", "a_minus", "a_diff"]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)

  skel_panel_ids = get_panel_ids(Λ)
  skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
  writevtk(Λ,dir*"/ambient_model_skeleton",cellfields=cellfields,append=false,geo_map=skel_geo_map)
end

################################################################################
## Advection tests
### check abs(v⋅n.plus) = abs(v⋅n.minus)
################################################################################
vecX(XYZ) = VectorValue(-XYZ[2],XYZ[3],0.0)
vX = panel_to_cartesian(tangent_vec(vecX))

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
U = TrialFESpace(V)

Ω_panel = Triangulation(panel_model)
panel_ids = get_panel_ids(panel_model)
Λ = SkeletonTriangulation(panel_model)
n_Λ = get_normal_vector(Λ)
pts = get_cell_points(Λ)

_vel = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
vel = interpolate(_vel,U)

diff_cf = (abs((vel⋅ n_Λ).minus) .- abs((vel⋅ n_Λ).plus))
@test all(all.(lazy_map(x->isless.(x, 1e-14), diff_cf(pts))))

if return_vtk
  skel_panel_ids = get_panel_ids(Λ)
  skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
  labels = ["upwind_plus","upwind_minus","upwind_diff"]
  panel_cfs = [abs((vel⋅ n_Λ).plus),abs((vel⋅ n_Λ).minus),abs((vel⋅ n_Λ).minus)-abs((vel⋅ n_Λ).plus)]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)
  writevtk(Λ,dir*"/ambient_model_skeleton", cellfields=cellfields,append=false,geo_map=skel_geo_map)
end


end
