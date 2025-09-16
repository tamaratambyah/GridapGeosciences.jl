panel_model = coarse_parametric_model()
panel_model = Gridap.Adaptivity.refine(panel_model)

panel_ids = get_panel_ids(panel_model)




Λ = SkeletonTriangulation(panel_model)

panelwise_cellfield(_analytic_inv_metric,Λ)

# _Λ = Λ.trian

# func = _analytic_inv_metric

function panelwise_cellfield(f::Function,atrian::AdaptedTriangulation)
  println("Adapted Boundary Triangulation")
  _cf = panelwise_cellfield(f,atrian.trian)
  if typeof(_cf)<:SkeletonPair
    println("skeleton ")
    plus = GenericCellField(get_data(_cf.plus),atrian,ReferenceDomain())
    minus = GenericCellField(get_data(_cf.minus),atrian,ReferenceDomain())
    return SkeletonPair(plus,minus)
  else
    return GenericCellField(get_data(_cf),atrian,ReferenceDomain())
  end

end

function panelwise_cellfield(f::Function,trian::BoundaryTriangulation)
  face_panel_ids = get_panel_ids(trian)

  @check length(face_panel_ids) == num_cells(trian) "\n Incorrect panel ids"

  cell_field = map(p->GenericField(f(p)),face_panel_ids)
  cf = CellData.GenericCellField(cell_field,trian,PhysicalDomain())

  _face_cf = boundary_cellfield(cf,trian)
  face_cf = GenericCellField(_face_cf,trian,ReferenceDomain())
  return face_cf
end

function boundary_cellfield(cf::CellField,trian::BoundaryTriangulation)
  glue = get_glue(trian,Val(2))

  panel_model = get_background_model(trian)
  cmap = get_cell_map(get_grid(panel_model))
  fcmap = lazy_map(Reindex(cmap),glue.tface_to_mface)

  ref_face_2_phys_cell_map = lazy_map(∘,fcmap,glue.tface_to_mface_map)
  face_cf = lazy_map(∘,get_data(cf),ref_face_2_phys_cell_map)
  return face_cf
end

function panelwise_cellfield(f::Function,trian::SkeletonTriangulation)
  ### plus
  panel_model = get_background_model(trian.plus)
  face_panel_ids = get_panel_ids(trian.plus)
  cell_field = map(p->GenericField(f(p)),face_panel_ids)
  cf_plus = CellData.GenericCellField(cell_field,trian.plus,PhysicalDomain())

  glue = get_glue(trian.plus,Val(2))

  cmap = get_cell_map(get_grid(panel_model))
  fcmap = lazy_map(Reindex(cmap),glue.tface_to_mface)

  ref_face_2_phys_cell_map = lazy_map(∘,fcmap,glue.tface_to_mface_map)
  _face_cf_plus = lazy_map(∘,get_data(cf_plus),ref_face_2_phys_cell_map)

  plus = GenericCellField(_face_cf_plus,trian,ReferenceDomain())



  #### minus
  face_panel_ids = get_panel_ids(trian.minus)
  cell_field = map(p->GenericField(f(p)),face_panel_ids)
  cf_minus = CellData.GenericCellField(cell_field,trian.minus,PhysicalDomain())

  glue = get_glue(trian.minus,Val(2))

  cmap = get_cell_map(get_grid(panel_model))
  fcmap = lazy_map(Reindex(cmap),glue.tface_to_mface)

  ref_face_2_phys_cell_map = lazy_map(∘,fcmap,glue.tface_to_mface_map)
  _face_cf_minus = lazy_map(∘,get_data(cf_minus),ref_face_2_phys_cell_map)

  minus = GenericCellField(_face_cf_minus,trian,ReferenceDomain())


  SkeletonPair(plus,minus)
end


### plus
trian = _Λ.plus
face_panel_ids = get_panel_ids(trian)
cf_plus = panelwise_cellfield(func,trian,face_panel_ids)

glue = get_glue(trian,Val(2))

cmap = get_cell_map(get_grid(panel_model))
fcmap = lazy_map(Reindex(cmap),glue.tface_to_mface)

ref_face_2_phys_cell_map = lazy_map(∘,fcmap,glue.tface_to_mface_map)
_face_cf_plus = lazy_map(∘,get_data(cf_plus),ref_face_2_phys_cell_map)

face_cf_plus = GenericCellField(_face_cf_plus,Λ,ReferenceDomain())

### minus
trian = _Λ.minus
face_panel_ids = get_panel_ids(trian)
cf_minus = panelwise_cellfield(func,trian,face_panel_ids)

glue = get_glue(trian,Val(2))

cmap = get_cell_map(get_grid(panel_model))
fcmap = lazy_map(Reindex(cmap),glue.tface_to_mface)

ref_face_2_phys_cell_map = lazy_map(∘,fcmap,glue.tface_to_mface_map)
_face_cf_minus = lazy_map(∘,get_data(cf_minus),ref_face_2_phys_cell_map)

face_cf_minus = GenericCellField(_face_cf_minus,Λ,ReferenceDomain())


# evaluate(_face_cf_plus[1],Point(1,))
# evaluate(_face_cf_minus[1],Point(1,))

skel_cf = SkeletonPair(face_cf_plus,face_cf_minus)
ginv_cf = SkeletonPair(face_cf_plus,face_cf_minus)


# pts = get_cell_points(Λ)
# skel_cf.plus(pts) .≈ skel_cf.minus(pts)
# skel_cf.plus(pts) .- skel_cf.minus(pts)

skel_ids = get_panel_ids(_Λ.plus)
skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_ids)
labels = ["plus","minus","plus-minus"]
panel_cfs = [skel_cf.plus,skel_cf.minus,skel_cf.plus-skel_cf.minus]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

writevtk(Λ,dir*"/ambient_model_skeleton", cellfields=cellfields,append=false,geo_map=skel_geo_map)
