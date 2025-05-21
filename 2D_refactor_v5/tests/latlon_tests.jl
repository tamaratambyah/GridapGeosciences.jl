using Gridap
include("../src/initialise.jl")


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
  ϕ
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


cdofs = get_cell_dof_values(uh)
for p in collect(1:6)
  panel_vals = cdofs[panel_ids.==p]
  println("---p = $p: Maximum = ", maximum(maximum.(panel_vals)) )
  println("---p = $p: Minimum = ", minimum(minimum.(panel_vals)) )
end

mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)
cf_mapped = lazy_map(Broadcasting(∘),get_data(uh),mapping)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)

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
    Jt_inv = evaluate(pinvJt(Jt),latlon_panelp)

    Jt_inv ⋅ uθϕ(latlon_panelp)
  end
end


################################################################################
##### Vector valued functions
################################################################################
function uθϕ_vector(θϕ)
  θ,ϕ = θϕ
  VectorValue(cos(θ),0.0)
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


writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)

# cdofs = get_cell_dof_values(uh)
# for p in collect(1:6)
#   panel_vals = cdofs[panel_ids.==p]
#   println("---p = $p: Maximum = ", maximum(maximum.(panel_vals)) )
#   println("---p = $p: Minimum = ", minimum(minimum.(panel_vals)) )
# end

mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

_uh = change_domain(uh,Ω_parametric,PhysicalDomain())
_cf_mapped = lazy_map(Broadcasting(push_∇),get_data(_uh),mapping)
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,inv_mapping)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)


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

# p1_bad_cells = get_faces(get_grid_topology(latlon_model),0,2)[bad_node]
# p1_bad_cells = [13, 14, 15, 16, 25, 26, 27, 28, 37, 38, 39, 40, 49, 50, 51, 52]
p1_bad_cells = [15, 16, 27, 28, 37, 38,  49, 50]
Nc = num_cells(latlon_model)
n_cells_per_panel = Nc/6
p2_bad_cells = Int(n_cells_per_panel) .+ p1_bad_cells
p5_bad_cells = Int(4*n_cells_per_panel) .+ p1_bad_cells

bad_cells = vcat(p2_bad_cells, p5_bad_cells)

good_cells = setdiff(collect(1:Nc),bad_cells)

Ω_view = Triangulation(ambient_model,good_cells)
cf_view = change_domain(cf_ambient,Ω_view,PhysicalDomain() )
writevtk(Ω_view,dir*"/ambient_view",cellfields=["u"=>cf_view],append=false)

################################################################################
##### Visualise at quad points
################################################################################
struct MyVisualizationGrid{Dc,Dp} <: Grid{Dc,Dp}
  sub_grid::UnstructuredGrid{Dc,Dp}
  sub_cell_to_cell::AbstractVector{<:Integer}
  cell_to_refpoints::AbstractVector{<:AbstractVector{<:Point}}
  cell_quad::AbstractVector{<:AbstractVector{<:Point}}
end

Gridap.Geometry.get_reffes(g::MyVisualizationGrid) = get_reffes(g.sub_grid)

Gridap.Geometry.get_cell_type(g::MyVisualizationGrid) = get_cell_type(g.sub_grid)

Gridap.Geometry.get_node_coordinates(g::MyVisualizationGrid) = get_node_coordinates(g.sub_grid)

Gridap.Geometry.get_cell_node_ids(g::MyVisualizationGrid) = get_cell_node_ids(g.sub_grid)

function MyVisualizationGrid(trian::Triangulation, ref_grids::AbstractArray{<:UnstructuredGrid})
  dΩ = Measure(trian,2)
  cell_quad = dΩ.quad.cell_point

  cell_to_ctype = collect1d(get_cell_type(trian))
  ctype_to_refpoints = map(get_node_coordinates, ref_grids)
  cell_to_refpoints = CompressedArray(ctype_to_refpoints,cell_to_ctype)
  cell_map = get_cell_map(trian)
  cell_to_points = lazy_map(evaluate,cell_map, cell_to_refpoints)

  node_to_coords, cell_to_offset = Gridap.Visualization._prepare_node_to_coords(cell_to_points)

  ctype_to_scell_to_snodes = map(get_cell_node_ids,ref_grids)

  sub_cell_to_nodes, sub_cell_to_cell = Gridap.Visualization._prepare_sub_cell_to_nodes(
    cell_to_ctype,ctype_to_scell_to_snodes,cell_to_offset)

  ctype_to_reffes = map(get_reffes,ref_grids)
  ctype_to_scell_type = map(get_cell_type,ref_grids)
  sctype_to_reffe, sub_cell_to_sctype = Gridap.Visualization._prepare_sctype_to_reffe(
    ctype_to_reffes,ctype_to_scell_type,cell_to_ctype)

  sub_grid = UnstructuredGrid(
    node_to_coords,
    sub_cell_to_nodes,
    sctype_to_reffe,
    sub_cell_to_sctype,
    NonOriented())

    MyVisualizationGrid(sub_grid,sub_cell_to_cell,cell_to_refpoints,cell_quad )

end


function Gridap.Visualization.visualization_data(
  trian::Triangulation, filebase::AbstractString;
  order=-1, nsubcells=-1, celldata=Dict(), cellfields=Dict())

  println("my vis")
  if order == -1 && nsubcells == -1
    # Use the given cells as visualization cells
    f = (reffe) -> UnstructuredGrid(reffe)
  elseif order != -1 && nsubcells == -1
    # Use cells of given order as visualization cells
    f = (reffe) -> UnstructuredGrid(LagrangianRefFE(Float64,get_polytope(reffe),order))
  elseif order == -1 && nsubcells != -1
    # Use linear sub-cells with nsubcells per direction
    f = (reffe) -> UnstructuredGrid(compute_reference_grid(reffe,nsubcells))
  else
    @unreachable "order and nsubcells kw-arguments can not be given at the same time"
  end

  ref_grids = map(f, get_reffes(trian))
  visgrid = MyVisualizationGrid(trian,ref_grids)

  cdata = Gridap.Visualization._prepare_cdata(celldata,visgrid.sub_cell_to_cell)
  pdata = Gridap.Visualization._prepare_pdata(trian,cellfields,visgrid.cell_quad)

  (Gridap.Visualization.VisualizationData(visgrid,filebase;celldata=cdata,nodaldata=pdata),)
end

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_view,dir*"/ambient_view",cellfields=["u"=>cf_view],append=false)
