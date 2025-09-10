
panel_model = coarse_parametric_model()
panel_model = Adaptivity.refine(panel_model)


panel_ids = get_panel_ids(panel_model)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

################################################################################
## Face normals
################################################################################
## Boundary trian
trian = BoundaryTriangulation(panel_model,1)
pts = get_cell_points(trian)

# get face normals in 3D
cell_vectors = Geometry.get_facet_normal(trian,cell_geo_map)
n_3D = get_normal_vector(trian,cell_vectors)

# push forward 2D chart normals to ambient space
n_mapped = pushforward_normal(trian)

# test equality
@test sum(n_mapped(pts) .≈ n_3D(pts)) == num_facets(panel_model)

################################################################################
### skeleton trian
trian = SkeletonTriangulation(panel_model)
pts = get_cell_points(trian)

# regular 2D normal in chart
n_2D = get_normal_vector(trian)

# get face normals in 3D
cell_vectors = Geometry.get_facet_normal(trian,cell_geo_map)
n_3D = get_normal_vector(trian,cell_vectors)

# push forward 2D chart normals to ambient space
n_mapped = pushforward_normal(trian)

# test equality
@test sum(n_mapped.plus(pts) .≈ n_3D.plus(pts)) == num_facets(panel_model)
@test sum(n_mapped.minus(pts) .≈ n_3D.minus(pts)) == num_facets(panel_model)

## plot normals on skeleton
panel_cfs = [n_3D.plus, n_3D.minus, n_3D.minus+n_3D.plus,
             n_2D.plus, n_2D.minus, n_2D.minus+n_2D.plus]
labels = ["amb_n_plus", "amb_n_minus", "amb_n_total",
          "chart_n_plus", "chart_n_minus", "chart_n_total"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

skel_panel_ids = get_panel_ids(trian)
skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
writevtk(trian,dir*"/ambient_model_skeleton",cellfields=cellfields,append=false,geo_map=skel_geo_map)

################################################################################
## Advection tests
################################################################################
### check sqrt(g) is continuous across skeleton
Λ = SkeletonTriangulation(panel_model)
skel_meas_cf = CellField(sqrtg,Λ)
pts = get_cell_points(Λ)

# test equality of plus and minus side
@test sum(skel_meas_cf.minus(pts) .≈ skel_meas_cf.plus(pts)) == num_facets(panel_model)

panel_cfs = [skel_meas_cf.plus, skel_meas_cf.minus, skel_meas_cf.minus-skel_meas_cf.plus]
labels = ["g_plus", "g_minus", "diff"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

skel_panel_ids = get_panel_ids(Λ)
skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
writevtk(Λ,dir*"/ambient_model_skeleton",cellfields=cellfields,append=false,geo_map=skel_geo_map)

################################################################################
### abs(v⋅n.plus) = abs(v⋅n.minus)
vecX(XYZ) = VectorValue(-XYZ[2],XYZ[3],0.0)
vX = panel_to_cartesian(tangent_vec(vecX))

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
U = TrialFESpace(V)

_vel = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
vel = interpolate(_vel,U)
labels = ["upwind_plus","upwind_minus","upwind_diff"]
panel_cfs = [abs((vel⋅ n_Λ).plus),abs((vel⋅ n_Λ).minus),abs((vel⋅ n_Λ).minus)-abs((vel⋅ n_Λ).plus)]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Λ,dir*"/ambient_model_skeleton", cellfields=cellfields,append=false,geo_map=skel_geo_map)
