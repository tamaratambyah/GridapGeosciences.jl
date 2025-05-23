using Gridap
include("../src/initialise.jl")
include("../src/Visualization/VisualizationData.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)
latlon_model = get_latlon_model(manifold_model)


Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
Ω_latlon = Triangulation(latlon_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)
pts_latlon = get_cell_points(Ω_latlon)

################################################################################
##### Map analytical scalar valued function from θ,ϕ -> α,β
##### 1. Map a point α,β -> θ,ϕ
##### 2. Evaluate g(α,β) = f(θ,ϕ)
################################################################################
function u_latlon(p::Int,uθϕ::Function)
  function _u(αβ)
    latlon_panel1 = evaluate(GnomonicMap(), αβ)
    sphere_panel1 = evaluate(Sigma(),latlon_panel1)
    sphere_panelp = evaluate(PanelRotationMap(r1p_3D[p]), sphere_panel1)
    latlon_panelp = evaluate(Sigma(),sphere_panelp)
    uθϕ(latlon_panelp)
  end
end



################################################################################
##### Scalar valued functions
################################################################################
function uθϕ_scalar(θϕ)
  θ,ϕ = θϕ
  rem2pi(θ,RoundNearest)
end

cell_field = map(p->GenericField(u_latlon(p,uθϕ_scalar)),panel_ids)
cf_parametric = CellData.GenericCellField(cell_field,Ω_parametric,PhysicalDomain())

dΩ = Measure(Ω_parametric,2)
H1 = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1), conformity=:H1)
biform(u,v) = ∫(u*v)dΩ
liform(v) = ∫(v*cf_parametric )dΩ

op = AffineFEOperator(biform,liform,H1,H1)
uh = solve(LUSolver(),op)

writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)

mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)
cf_mapped = lazy_map(Broadcasting(∘),get_data(uh),mapping)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )

writevtk(cubedsphere,Ω_ambient,dir*"/ambient_quad",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)

### compare values
cdofs = get_cell_dof_values(uh)
for p in collect(1:6)
  panel_vals = cdofs[panel_ids.==p]
  println("---p = $p: Maximum = ", maximum(maximum.(panel_vals)) )
  println("---p = $p: Minimum = ", minimum(minimum.(panel_vals)) )
end

ambient_vals = cf_ambient(pts_ambient)
for p in collect(1:6)
  panel_vals = ambient_vals[panel_ids.==p]
  println("---p = $p: Maximum = ", maximum(maximum.(panel_vals)) )
  println("---p = $p: Minimum = ", minimum(minimum.(panel_vals)) )
end

################################################################################
##### Map analytical vector valued function from θ,ϕ -> α,β
##### 1. Map a point α,β -> θ,ϕ.
##### 2. Compute u(θ,ϕ)
##### 2. Compute J cooresponding to θ,ϕ -> α,β
##### 2. Evaluate v(α,β) = J^{-T} ⋅ u(θ,ϕ)
################################################################################

function u_latlon_vector(p::Int,uθϕ::Function)
  function _u(αβ)
    latlon_panel1 = evaluate(GnomonicMap(), αβ)
    sphere_panel1 = evaluate(Sigma(),latlon_panel1)
    sphere_panelp = evaluate(PanelRotationMap(r1p_3D[p]), sphere_panel1)
    latlon_panelp = evaluate(Sigma(),sphere_panelp)

    # cmap = InvSigmaField(r) ∘ PanelRotationField(r1p_3D[p]) ∘ SigmaField(r) ∘ GnomonicField()
    cmap = InvGnomonicField() ∘ InvSigmaField(r)  ∘ PanelRotationField(rp1_3D[p]) ∘ SigmaField(r)
    Jt = ∇(cmap)
    det_J = meas(Jt)
    Jt_inv = evaluate(pinvJt(Jt),latlon_panelp)

    det_J * Jt_inv ⋅ uθϕ(latlon_panelp)
  end
end


################################################################################
##### Vector valued functions
################################################################################
function uθϕ_vector(θϕ)
  θ,ϕ = θϕ
  VectorValue(cos(θ),0.0)
  VectorValue(X,Y,Z)
end



dΩ = Measure(Ω_parametric,2)
RT = FESpace(manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)

## interpolate everywhere
free_values = zero_free_values(RT)
dirichlet_values = zero_dirichlet_values(RT)

## get cell vals
s = get_fe_dof_basis(RT)
trian = get_triangulation(s)
cell_field = map(p->GenericField(u_latlon_vector(p,uθϕ_vector)),panel_ids)
f = CellData.GenericCellField(cell_field,trian,PhysicalDomain())
cell_vals = s(f)

## interpolate!
gather_free_and_dirichlet_values!(free_values,dirichlet_values,RT,cell_vals)
uh = FEFunction(RT,free_values,dirichlet_values)


writevtk(Ω_parametcubedsphere,ric,dir*"/parametric",cellfields=["u"=>uh],append=false)

# cdofs = get_cell_dof_values(uh)
# for p in collect(1:6)
#   panel_vals = cdofs[panel_ids.==p]
#   println("---p = $p: Maximum = ", maximum(maximum.(panel_vals)) )
#   println("---p = $p: Minimum = ", minimum(minimum.(panel_vals)) )
# end

mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

Jt = lazy_map(Broadcasting(gradient),mapping)
Jt_inv = lazy_map(Broadcasting(pinvJt),Jt)
det_J = lazy_map(Operation(Mymeas),Jt_inv)
_det_J = lazy_map(Broadcasting(∘),det_J,inv_mapping)

cf_det_J = CellData.GenericCellField(_det_J,Ω_ambient,PhysicalDomain() )
writevtk(cubedsphere,Ω_ambient,dir*"/ambient_vector_det",cellfields=["u"=>cf_det_J],append=false)

function u_func(X)
  VectorValue(X[1],X[2],X[3])
end
writevtk(Ω_ambient,dir*"/ambient_vector_det",cellfields=["u"=>u_func],append=false)


function f(XY)

  if pole
    VectorValue(XY[1],0.0,0.0)
  else
    VectorValue(XY[1],XY[2],XY[3])
  end


end















function Mymeas(Jt::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  J = transpose(Jt)
  (det(Jt⋅J))
end


change = lazy_map(*,det_J,Jt_inv)

_uh = change_domain(uh,Ω_parametric,PhysicalDomain())
# _cf_mapped = lazy_map(Broadcasting(push_∇),get_data(_uh),mapping)
_cf_mapped = lazy_map(Broadcasting(push_∇),get_data(_uh),mapping)
_cf_mapped2 = lazy_map(*,det_J,_cf_mapped)
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped2,inv_mapping)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)
writevtk(cubedsphere,Ω_ambient,dir*"/ambient_vector",cellfields=["u"=>cf_ambient, "f"=>f_cf_ambient],append=false)


f_cf_mapped = lazy_map(Broadcasting(push_∇),get_data(f),mapping)
f_cf_mapped = lazy_map(Broadcasting(∘),f_cf_mapped,inv_mapping)
f_cf_ambient = CellData.GenericCellField(f_cf_mapped,Ω_ambient,PhysicalDomain() )
################################################################################
##### Visualisation hack: remove the bad cells
################################################################################
#### find bad cells
latlon_coords = get_cell_coordinates(latlon_model)
# find the zero on panel 1
p1 = latlon_coords[panel_ids.==1]

found = false
cell = 0
idx = 0
for i in 1:length(p1)
  pp = p1[i]
  for j in 1:length(pp)
    if pp[j] == VectorValue(0.0,0.0)
      println(j)
      cell = i
      idx = j
      found = true
      break
    end
  end
  if found
    break
  end
end

bad_node = get_cell_node_ids(latlon_model)[cell][idx]

p1_bad_cells = get_faces(get_grid_topology(latlon_model),0,2)[bad_node]
# p1_bad_cells = [13, 14, 15, 16, 25, 26, 27, 28, 37, 38, 39, 40, 49, 50, 51, 52]
# p1_bad_cells = [15, 16, 27, 28, 37, 38,  49, 50]
Nc = num_cells(latlon_model)
n_cells_per_panel = Nc/6
p2_bad_cells = Int(n_cells_per_panel) .+ p1_bad_cells
p5_bad_cells = Int(4*n_cells_per_panel) .+ p1_bad_cells

bad_cells = vcat(p2_bad_cells, p5_bad_cells)

good_cells = setdiff(collect(1:Nc),bad_cells)

Ω_view = Triangulation(ambient_model,good_cells)
cf_view = change_domain(cf_ambient,Ω_view,PhysicalDomain() )
f_cf_view = change_domain(f_cf_ambient,Ω_view,PhysicalDomain() )
writevtk(Ω_view,dir*"/ambient_view2",cellfields=["u"=>cf_view, "f"=>f_cf_view],append=false,nsubcells=8)



writevtk(Ω_ambient,dir*"/ambient_",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_view,dir*"/ambient_view",cellfields=["u"=>cf_view],append=false)
