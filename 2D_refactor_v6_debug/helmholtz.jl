include("../2D_refactor_v6/overloads.jl")

degree = 10


dΩ = Measure(Ω_panel,degree)


cell_field = map(p->GenericField(measure(p)),panel_ids)
meas_cf =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

cell_field = map(p->GenericField(inv_metric(p)),panel_ids)
inv_metric_cf = CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,2); conformity=:H1)
U = TrialFESpace(V)

uh_sin = interpolate(usin_cf_panel,U)

u = -uh_sin
slap = 2*uh_sin

rhs = u + slap
sum(∫( rhs*meas_cf  )dΩ)


poisson_biform(u,v) = ∫(u*v*meas_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
poisson_liform(v) = ∫(  (rhs*v)*meas_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-u,dΩ)

for pid in collect(1:6)
  mask = panel_ids.==pid
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/u$pid",cellfields=["u"=>u,"uh"=>uh,"e"=>u-uh],append=false)
end


### map to sphere
inv_f = lazy_map(p->InverseMap(p),panel_ids)

_uh = change_domain(uh,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh),inv_f)
uh_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field

_u = change_domain(u,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_u),inv_f)
u_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field

writevtk(Triangulation(ambient_model),dir*"/ambient_model_helmholtz",
          cellfields=["u"=>u_ambient,"uh"=>uh_ambient,"e"=>u_ambient-uh_ambient
          ],append=false)



########### possion
V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,1); conformity=:H1, constraint=:zeromean)
U = TrialFESpace(V)

uh_sin = interpolate(usin_cf_panel,U)

u = -uh_sin
slap = 2*uh_sin
rhs = slap

sum( ∫(  rhs)dΩ )

poisson_biform(u,v) = ∫( ((basis_vectors ⋅ inv_metric_cf ⋅ gradient(v)) ⋅ (basis_vectors ⋅ inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
poisson_liform(v) = ∫(  (rhs*v)*meas_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

A_mat = get_matrix(op)
using LinearAlgebra
eigvals(Array(A_mat))


Λ = ConstantFESpace(panel_model)
M = TrialFESpace(Λ)
X = MultiFieldFESpace([U,M])
Y = MultiFieldFESpace([V,Λ])
poisson_biformX((u,μ),(v,λ)) = poisson_biform(u,v)  + ∫( (v*μ)*meas_cf )dΩ + ∫((λ*u)*meas_cf)dΩ
poisson_liformY((v,λ)) =   poisson_liform(v)  + ∫( (λ*uh_sin)*meas_cf )dΩ

op = AffineFEOperator(poisson_biformX,poisson_liformY,X,Y)
uhX,μh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-u,dΩ)


for p in collect(1:6)
  mask = panel_ids.==p
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/panel$(p)_model",
          cellfields=["uh_sin"=>uh_sin, "uh"=>uhX
          ],append=false)
end
