## make models
panel_model = coarse_parametric_model()
panel_model = Gridap.Adaptivity.refine(panel_model)
panel_model = Gridap.Adaptivity.refine(panel_model)



## arbitary functions
# depth(XYZ) = 1.0 + 0.1*exp(-( XYZ[2]^2 + XYZ[3]^2 ) )
# velocity(XYZ) = VectorValue(-XYZ[2],XYZ[1],0.0)
# coriolis(XYZ) = 1e-5
# h = panel_to_cartesian(depth)
# vecX = velocity
# vX = panel_to_cartesian(tangent_vec(vecX))
# f = panel_to_cartesian(coriolis)

## williamson functions
ζ,u0,ω = 0.0, 0.1, 1e-5
h = panel_to_latlon(hWilliamson(ζ,u0,ω))
vecX = vec_cartesian_to_latlon(vWilliamson(ζ,u0,ω))
vX = panel_to_cartesian(tangent_vec(vecX))
f = panel_to_latlon(fWilliamson(ζ,u0,ω))

gravity = 1.0

p_fe = 1


panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*(p_fe+1))

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)


R = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
H = TrialFESpace(R)


Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])


covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)

u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)




# weak forms
detg_cf = CellField(detg,Ω_panel)
metric_cf = CellField(analytic_metric,Ω_panel)
inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)
grad_meas_cf = CellField(grad_meas,Ω_panel)


## initial conditions
u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_contra_h = interpolate(u_contra_cf,U)
u_proj_h = covarient_basis_cf ⋅ u_contra_h

h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
h_h = interpolate(h_cf,P)


# mass flux
biformF(F,v) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ
liformF(v) = ∫( h_h*(u_contra_h⋅(metric_cf⋅v))*meas_cf   )dΩ
op = AffineFEOperator(biformF,liformF,U,V)
Fh = solve(LUSolver(),op)

### equation for depth:
rhs_h = h_h + 1/meas_cf*( Fh⋅grad_meas_cf + meas_cf*(∇⋅Fh)   )
biform_p(p,r) = ∫( (p*r)*meas_cf )dΩ
liform_p(r) = ( ∫( (rhs_h*r)*meas_cf )dΩ
              - ∫( r*(Fh⋅grad_meas_cf + meas_cf*(∇⋅Fh) )  )dΩ
              )
op = AffineFEOperator(biform_p,liform_p,P,Q)
ph = solve(LUSolver(),op)
e_p = l2((h_cf - ph)*meas_cf,dΩ) # error in depth

panel_cfs = [ph, h_cf,ph-h_cf]
labels = ["ph","p","ep"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)





# Bernoulli potential
biformΦ(Φ,r) = ∫( Φ*r*meas_cf  )dΩ
liformΦ(r) = ∫( gravity*h_h*r*meas_cf  )dΩ + ∫( 0.5*( u_contra_h ⋅(metric_cf⋅u_contra_h) )r*meas_cf  )dΩ
op = AffineFEOperator(biformΦ,liformΦ,P,Q)
Φh = solve(LUSolver(),op)



# vorticity
_analytic_perp_matrix(αβ) = TensorValue{2,2}( -F(αβ), E(αβ), -G(αβ), F(αβ) )
perp_matrix_cf = CellField(_analytic_perp_matrix,Ω_panel)

biformq(q,r) = ∫( q*h_h*r*meas_cf  )dΩ
liformq(r) = ∫( cor_cf*r*meas_cf  )dΩ + ∫( (perp_matrix_cf⋅u_contra_h)⋅∇(r)  )dΩ
op = AffineFEOperator(biformq,liformq,H,R)
qh = solve(LUSolver(),op)


Fperph = 1/meas_cf*(perp_matrix_cf⋅Fh)
rhs_u = u_contra_h + qh*Fperph + (  inv_metric_cf⋅gradient(Φh) )


#### solve velocity
biform_u(u,v) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ
liform_u(v) = ( ∫( rhs_u⋅(metric_cf⋅v)*meas_cf )dΩ
                + ∫( Φh*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                 - ∫( qh*( (perp_matrix_cf⋅Fh) ⋅(metric_cf ⋅v))   )dΩ
                  )
op = AffineFEOperator(biform_u,liform_u,U,V)
uh = solve(LUSolver(),op)

uh_proj = covarient_basis_cf ⋅ uh

e_u = l2( ( uh-u_contra_h  )*meas_cf,dΩ  )
e_u = l2( (u_proj_h - uh_proj)*meas_cf,dΩ) # error in physical velocity u


panel_cfs = [uh_proj, u_proj_h, uh_proj-u_proj_h]
labels = ["uh","u","eu"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)




## mass matrices
biform_u((u,p),(v,r)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ
biform_p((u,p),(v,r)) = ∫( (p*r)*meas_cf )dΩ
biformX((u,p),(v,r)) = biform_u((u,p),(v,r)) + biform_p((u,p),(v,r))
A = assemble_matrix(biformX,X,Y)

## rhs vectors
liform1((v,r)) = ∫( rhs_u⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( (rhs_h*r)*meas_cf )dΩ
liform2((v,r)) =  ( ∫( Φh*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                  - ∫( qh*( Fperph ⋅(metric_cf ⋅v)) * meas_cf  )dΩ
                  - ∫( r*(Fh⋅grad_meas_cf + meas_cf*(∇⋅Fh) )  )dΩ
                  )
liformX((v,r)) = liform1((v,r)) + liform2((v,r))

b1 = assemble_vector(liform1,Y)
b = assemble_vector(liform2,Y) + b1

ns = Gridap.Algebra.LUNumericalSetup(lu(A))
x = allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)
xh = FEFunction(X,x)
uh,ph = xh

# op = AffineFEOperator(biformX,liformX,X,Y)
# uh,ph = solve(LUSolver(),op)

uh_proj = covarient_basis_cf ⋅ uh

e_u = l2( ( uh-u_contra_h  )*meas_cf,dΩ  )
e_u = l2( (u_proj_h - uh_proj)*meas_cf,dΩ) # error in physical velocity u
e_p = l2((h_cf - ph)*meas_cf,dΩ) # error in depth


cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
panel_cfs = [ph, uh_proj, uh_proj-u_proj_cf,ph-h_cf]
labels = ["p","u_proj","eu","ep","rhs_s"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)
