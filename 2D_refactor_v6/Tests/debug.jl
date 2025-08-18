
####################### refine debug
panel_model = coarse_parametric_model()

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), get_panel_ids(panel_model))
writevtk(Triangulation(panel_model),dir*"/ambient_grid",append=false,geo_map=cell_geo_map)

amodel = Gridap.Adaptivity.refine(panel_model)
r_panel_ids = get_panel_ids(amodel)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), r_panel_ids)
writevtk(Triangulation(amodel),dir*"/ref_ambient_grid",append=false,geo_map=cell_geo_map)
panel_model = amodel
