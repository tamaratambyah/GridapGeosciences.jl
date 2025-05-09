using Gridap
include("../src/initialise.jl")

#############
a = π/4
u_parametric(x) = (a+x[1])*(a-x[1]) + (a+x[2])*(a-x[2])
p = 2

######### cubed sphere model
coarse_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
model = Adaptivity.refine(coarse_model)
ref_model = Adaptivity.refine(Adaptivity.refine(model))

manifold_model = ref_model
panel_ids = get_panel_ids(manifold_model)

############# Surface metric and quadrature
order = 2*(p+1)
Ω_parametric = Triangulation(manifold_model)

m = Metric(cubedsphere,Ω_parametric)
dΩg = Measure(m,Ω_parametric,order)

############# FE problem
u_cf = CellField(u_parametric,Ω_parametric)
f_cf = -1.0*surface_laplacian(u_cf,m)

V = FESpace(manifold_model, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

biform(u,v) = ∫(  surface_gradient(v,m) ⋅ surface_gradient(u,m)  )dΩg
liform(v)   = ∫( v*f_cf )dΩg

op = AffineFEOperator(biform,liform,U,V)
ls = LUSolver()
uh = solve(ls,op)

e = uh-u_cf
l2(e,dΩg)

A = get_matrix(op)
b = get_vector(op)




# plot in parametric space
writevtk(Ω_parametric,dir*"/possion_parametric",
        cellfields=["uh"=>uh,"u"=>u_cf, "e"=>e, "f"=>f_cf],append=false)
