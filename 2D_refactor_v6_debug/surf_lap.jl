function rho(αβ)
  α,β = αβ
  sqrt(1 + (tan(α))^2 + (tan(β))^2 )
end

function rho3(αβ)
  α,β = αβ
  ( sqrt(1 + (tan(α))^2 + (tan(β))^2 ) )^3
end

function dXda(αβ)
  α,β = αβ
  - tan(α)*(sec(α))^2 / ( rho3(αβ)  )
end

function dXdb(αβ)
  α,β = αβ
  - tan(β)*(sec(β))^2 /(  rho3(αβ)  )
end

function dYda(αβ)
  α,β = αβ
  dXda(αβ)*tan(α) + 1/rho(αβ)*(sec(α))^2
end

function dYdb(αβ)
  α,β = αβ
  dXdb(αβ)*tan(α)
end

function dZda(αβ)
  α,β = αβ
  dXda(αβ)*tan(β)
end

function dZdb(αβ)
  α,β = αβ
  dXdb(αβ)*tan(β) + 1/rho(αβ)*(sec(β))^2
end

E(αβ) = dXda(αβ)*dXda(αβ) + dYda(αβ)*dYda(αβ) + dZda(αβ)*dZda(αβ)
F(αβ) = dXda(αβ)*dXdb(αβ) + dYda(αβ)*dYdb(αβ) + dZda(αβ)*dZdb(αβ)
G(αβ) = dXdb(αβ)*dXdb(αβ) + dYdb(αβ)*dYdb(αβ) + dZdb(αβ)*dZdb(αβ)

detg(αβ) = E(αβ)*G(αβ) - F(αβ)*F(αβ)
sqrtg(αβ) = sqrt( E(αβ)*G(αβ) - F(αβ)*F(αβ) )
analytic_metric(αβ) = TensorValue{2,2}(E(αβ),F(αβ),F(αβ),G(αβ))
analytic_inv_metric(αβ) =  TensorValue{2,2}(G(αβ)/detg(αβ),-F(αβ)/detg(αβ),-F(αβ)/detg(αβ),E(αβ)/detg(αβ))

## f = XYZ
# function f(p)
#   function _f(αβ)
#     α,β = αβ
#     if p == 1 || p == 5 || p == 6
#       return RADIUS^3/rho3(αβ)*tan(α)*tan(β)
#     else
#       return -RADIUS^3/rho3(αβ)*tan(α)*tan(β)
#     end
#   end
# end

### f = sin(ϕ) = Z
# function f(p)
#   function _f(αβ)
#     α,β = αβ
#     if p == 1
#       return 1/rho(αβ)*tan(β)
#     elseif p == 2
#       return 1/rho(αβ)
#     elseif p == 3
#       return 1/rho(αβ)*tan(β)
#     elseif p == 4
#       return 1/rho(αβ)*tan(α)
#     elseif p == 5
#       return -1/rho(αβ)
#     elseif p == 6
#       return 1/rho(αβ)*tan(α)
#     end
#   end
# end

### f = cos(ϕ) = cos(Z)
function f(p)
  function _f(αβ)
    α,β = αβ
    if p == 1
      return cos(1/rho(αβ)*tan(β))
    elseif p == 2
      return cos(1/rho(αβ))
    elseif p == 3
      return cos(1/rho(αβ)*tan(β))
    elseif p == 4
      return cos(1/rho(αβ)*tan(α))
    elseif p == 5
      return cos(-1/rho(αβ))
    elseif p == 6
      return cos(1/rho(αβ)*tan(α))
    end
  end
end

dfda(p) = αβ -> (gradient(f(p))(αβ))[1]
dfdb(p) = αβ -> (gradient(f(p))(αβ))[2]

gradient(f(1))(Point(1,1))

surflap(1)(Point(1,1))


w1(p) = αβ -> 1/sqrtg(αβ) * G(αβ)*dfda(p)(αβ) - 1/sqrtg(αβ)*F(αβ)*dfdb(p)(αβ)
w2(p) = αβ -> -1/sqrtg(αβ) * F(αβ)*dfda(p)(αβ) + 1/sqrtg(αβ)*E(αβ)*dfdb(p)(αβ)

w(p) = αβ -> VectorValue(w1(p)(αβ),w2(p)(αβ))

surflap(p) = αβ -> 1/sqrtg(αβ)*(divergence(w(p))(αβ))

rhs_func(p) = αβ -> f(p)(αβ) + surflap(p)(αβ)

Ω_panel = Triangulation(panel_model)
Ω_sphere = Triangulation(ambient_model)


cell_field = map(p->GenericField(f(p)),panel_ids)
f_cf =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

cell_field = map(p->GenericField(surflap(p)),panel_ids)
surflap_cf =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

cell_field = map(p->GenericField(rhs_func(p)),panel_ids)
rhs_cf =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())


inv_f = lazy_map(p->InverseMap(p),panel_ids)

cf_mapped = lazy_map(Broadcasting(∘),get_data(surflap_cf),inv_f)
slap_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field

cf_mapped = lazy_map(Broadcasting(∘),get_data(f_cf),inv_f)
f_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field

cf_mapped = lazy_map(Broadcasting(∘),get_data(rhs_cf),inv_f)
rhs_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field


writevtk(Triangulation(ambient_model),dir*"/ambient_model_slap",
          cellfields=["slap"=>slap_ambient,"f"=>f_ambient, "rhs"=>rhs_ambient
          ],append=false)




### helmholtz
dΩ = Measure(Ω_panel,10)


meas_cf = CellField(sqrtg,Ω_panel)
inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)

V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,2); conformity=:H1)
U = TrialFESpace(V)


sum(∫( rhs_cf*meas_cf  )dΩ)


poisson_biform(u,v) = ∫(u*v*meas_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-f_cf,dΩ)


_uh = change_domain(uh,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh),inv_f)
uh_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # am

writevtk(Triangulation(ambient_model),dir*"/ambient_model_slap",
          cellfields=["slap"=>slap_ambient,"f"=>f_ambient, "rhs"=>rhs_ambient,
          "uh"=>uh_ambient, "e"=>uh_ambient-f_ambient
          ],append=false)
