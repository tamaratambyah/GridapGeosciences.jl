include("../2D_refactor_v6/overloads.jl")
p_fe = 2
degree = 6*p_fe

Ω = Triangulation(panel_model)
dΩ = Measure(Ω,degree)

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
U = TrialFESpace(V)

uX(XYZ) =  XYZ[1]*XYZ[2]*XYZ[3]
cell_field = map(p->GenericField(uX_panel(p,uX)),panel_ids)
ucf =  CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

cell_field = map(p->GenericField(measure(p)),panel_ids)
meas_cf =  CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

cell_field = map(p->GenericField(inv_metric(p)),panel_ids)
inv_metric_cf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

uh = interpolate(ucf,U)

function _uX_panel(αβ,p)
  XYZ = forward_map(αβ,p)
  uX(XYZ)
end

function slap(p)
  function _slap(αβ)
    gu(αβ) = gradient(_uX_panel(αβ,p))
      inv_metric(αβ,p) ⋅ gu(αβ)
  end
end

cell_field = map(p->GenericField(slap(p)),panel_ids)
surflap = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())
surflap(get_cell_points(Ω))

gradient(uh)
laplacian(uh)
divergence(  gradient(uh)  )

rhs = ucf + surflap
sum(∫( rhs*meas_cf  )dΩ)





l2_biform(u,v) = ∫(u*v*meas_cf)dΩ
l2_liform(v) = ∫(  (ucf*v)*meas_cf )dΩ
op = AffineFEOperator(l2_biform,l2_liform,U,V)
uh_l2 = solve(LUSolver(),op)

poisson_biform(u,v) = ∫(u*v*meas_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
poisson_liform(v) = ∫(  (rhs*v)*meas_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-ucf,dΩ)

for pid in collect(1:6)
  mask = panel_ids.==pid
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/u$pid",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh],append=false)
end
