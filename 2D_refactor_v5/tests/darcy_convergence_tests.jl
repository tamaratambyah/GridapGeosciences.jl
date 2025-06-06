using Gridap
include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

###################### FE SOLVE
p = 1
degree = 2*p+1

# g_ambient(x) = sin(0.5*π*x[2]) # x[1]*x[2]*x[3]
# scalar_cf = map(p->GenericField(u_scalar_ambient2parametric(p,g_ambient)),panel_ids)
# g = CellData.GenericCellField(scalar_cf,Ω_parametric,PhysicalDomain())


RT = ReferenceFE(raviart_thomas,Float64,p)
DG = ReferenceFE(lagrangian,Float64,p)

V = FESpace(manifold_model, RT, conformity=:Hdiv)
Q = FESpace(manifold_model, DG, conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U, P])
Y = MultiFieldFESpace([V, Q])


m = Metric(cubedsphere,Ω_parametric)
dΩg = Measure(m,Ω_parametric,degree)


# darcy_biform((u,p),(v,q)) = ∫( v⋅u + p*q - p*(surface_divergence(v,m)) )dΩg# + q*(surface_divergence(u,m)) )dΩg
darcy_biform((u,p),(v,q)) = ∫( (p*( (1/m.sq_meas * divergence(m.sq_meas * v) ) ))*m.sq_meas )dΩ_parametric# + q*(surface_divergence(u,m)) )dΩg
darcy_liform((v,q)) = ∫(  q*0.0 )dΩg



# op = AffineFEOperator(darcy_biform,darcy_liform,X,Y)

### low levels
assem = SparseMatrixAssembler(X,Y)
dx = get_trial_fe_basis(X)
dy = get_fe_basis(Y)

matdata = collect_cell_matrix(X, Y, darcy_biform(dx,dy))
vecdata = collect_cell_vector(Y,darcy_liform(dy))

b = assemble_vector(assem,vecdata)
A = assemble_matrix(assem,matdata)

### low level matrix assembly
w = []
_r = []
c = []
dc = darcy_biform(dx,dy)
  for strian in get_domains(a)
    scell_mat = get_contribution(a,strian)
    cell_mat, trian = move_contributions(scell_mat,strian)
    @assert ndims(eltype(cell_mat)) == 2
    cell_mat_c = attach_constraints_cols(trial,cell_mat,trian)
    cell_mat_rc = attach_constraints_rows(test,cell_mat_c,trian)
    rows = get_cell_dof_ids(test,trian)
    cols = get_cell_dof_ids(trial,trian)
    #push!(w,compress_contributions(cell_mat_rc,trian))
    #push!(r,compress_ids(rows,trian))
    #push!(c,compress_ids(cols,trian))
    push!(w,cell_mat_rc)
    push!(_r,rows)
    push!(c,cols)
  end
  (w,_r,c)



xh = solve(LUSolver(),op)
uh, ph = xh

uh_ambient = parametric_cf_2_ambient_vector(manifold_model,uh)
ss, ph_ambient = parametric_cf_2_ambient(manifold_model,degree,ph)

writevtk(Ω_ambient,dir*"/darcy_ambient",
        cellfields=["uh"=>uh_ambient,"ph"=>ph_ambient],append=false)

# e = g-rh
# sqrt(sum(∫(e*e)dΩg))
