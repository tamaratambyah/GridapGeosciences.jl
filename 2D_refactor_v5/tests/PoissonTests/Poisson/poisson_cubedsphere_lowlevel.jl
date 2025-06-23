using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../../../src/initialise.jl")
include("../poisson_helpers.jl")



#### Analytic solution with zero mean
p = 2
degree = 20#2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

n = ns[3]

global RADIUS = 1.0*sqrt(3.0)


################################################################################
#### Analytic parametric space for a single panel
################################################################################

model = CartesianDiscreteModel((-π/4,π/4, -π/4,π/4), (n,n))
Ω = Triangulation(model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

function uex(x)
  if x[2] < 0.0
    return -x[2]*(x[2] + π/4)
  else
    return x[2]*(x[2] - π/4)
  end

end

# uex(x) = (x[2]+π/4)*(x[2] - π/4)*(x[1]+π/4)*(x[1]-π/4)
writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex],append=false)

ucf = CellField(uex,Ω)

# # check zero mean
# sum(∫(ucf)dΩ  )
# sum(∫(ucf)dΩg  )

# # check compatibility
# sum(∫( surface_laplacian(ucf,m))dΩ  )
# sum(∫( surface_laplacian(ucf,m))dΩg  )


# ## zero mean in FE space
# _rhs = ucf + -1.0*surface_laplacian(ucf,m)
# sum(∫(_rhs )dΩg)
# sum(∫(_rhs )dΩ)

# V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
# U = TrialFESpace(V)

# Λ = ConstantFESpace(model)
# M = TrialFESpace(Λ)

# X = MultiFieldFESpace([U,M])
# Y = MultiFieldFESpace([V,Λ])



# poisson_biform(u,v) = ∫(u*v)dΩg + ∫( surface_gradient(u,m)⋅gradient(v) )dΩg
# poisson_liform(v) =  ∫( (_rhs*v) )dΩg

# op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
# A = get_matrix(op)
# b = get_vector(op)
# eigvals(Array(A))
# _vec = A*ones(size(b))
# norm(A*ones(size(b)))
# sum(b)

# uh = solve(LUSolver(),op)


# e = l2(uh-ucf,dΩ)
# eg = l2(uh-ucf,dΩg)

# # poisson_biform((u,μ),(v,λ)) = ∫( u*v )dΩg + ∫( surface_gradient(u,m)⋅gradient(v) )dΩg   + ∫(v*μ)dΩg + ∫(λ*u)dΩg
# # poisson_liform((v,λ)) =  ∫( (_rhs*v) )dΩg
# # op = AffineFEOperator(poisson_biform,poisson_liform,X,Y)
# # uh,μh = solve(LUSolver(),op)



# b = get_vector(op)
# A = get_matrix(op)
# using LinearAlgebra
# eigvals(Array(A))
# _vec = A*ones(size(b))
# norm(A*ones(size(b)))
# sum(b)


# # eigvecs(Array(A))

# # labels = get_face_labeling(model)
# # add_tag_from_tags!(labels,"bb",[7])

# # Γ = BoundaryTriangulation(model,tags=["boundary"])
# # get_cell_dof_ids(V,Γ)

# # sum(b)
# # println("Compatibility: ", sum(b))


# e = l2(uh-ucf,dΩ)
# eg = l2(uh-ucf,dΩg)


# writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex,"uh"=>uh,"e"=>uex-uh],append=false)

# 1;
# # writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex],append=false)

p=4
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2)
U = TrialFESpace(V)

T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
S = TrialFESpace(T)

Λ = ConstantFESpace(model)
M = TrialFESpace(Λ)

X = MultiFieldFESpace([S,U,M])
Y = MultiFieldFESpace([T,V,Λ])

_X = MultiFieldFESpace([S,M])
_Y = MultiFieldFESpace([T,Λ])

#### With metric: Method 4 - mixed form -- with lagrange multiplers
### Sigma_exact
_sigma_ex = surface_gradient(ucf,m)

biformS((s,μ),(t,λ)) = ∫( s⋅t )dΩg + ∫( surface_divergence(s,m)*λ )dΩg  + ∫( surface_divergence(t,m)*μ )dΩg
liformS((t,λ)) = ∫( _sigma_ex ⋅ t )dΩg
op = AffineFEOperator(biformS,liformS,_X,_Y)
sigma_exh,μh = solve(LUSolver(),op)

# biformS(s,t) = ∫( s⋅t )dΩg
# liformS(t) = ∫( _sigma_ex ⋅ t )dΩg
# op = AffineFEOperator(biformS,liformS,S,T)
# sigma_exh = solve(LUSolver(),op)

# writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex,"s"=>sigma_ex,"sh"=>sigma_exh],append=false)


### dual form
# _rhs = -1.0*surface_divergence(sigma_exh,m)
_rhs = -1.0*surface_divergence(_sigma_ex,m)

biformX((s,u,μ),(t,v,λ)) = (  ∫( s⋅t  )dΩg + ∫( wave_divergence(t,m)*u )dΩ
                            + ∫( surface_divergence(s,m)*v  )dΩg
                            + ∫(v*μ)dΩg + ∫(λ*u)dΩg
                      )
liformY((t,v,λ)) = ∫( -(_rhs*v) )dΩg  + ∫(λ*ucf)dΩg

op = AffineFEOperator(biformX,liformY,X,Y)
sh,uh,μh = solve(LUSolver(),op)

#### Compute errors
e = sum(∫((uh-ucf)⊙(uh-ucf))dΩ)
eg = sum(∫((uh-ucf)⊙(uh-ucf))dΩg)

writevtk(Ω,dir*"/poisson",cellfields=["u"=>ucf,"uh"=>uh,"eu"=>ucf-uh,
   ],append=false)





################################################################################
#### With cubed sphere mesh
################################################################################

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

Ω = Triangulation(manifold_model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

function uex_p(p)
  # println("Panel: $p")
  function _u(αβ)
    if αβ[1] < 0.0
      # println("left ",-αβ[1]*(αβ[1] + π/4))
      return -αβ[1]*(αβ[1] + π/4)
    else
      # println("right", αβ[1]*(αβ[1] - π/4))
      return αβ[1]*(αβ[1] - π/4)
    end
  end
end

# _panel_ids = panel_ids[panel_ids.==1]

cell_field = map(p->GenericField(uex_p(p)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

pts = get_cell_points(Ω )
ucf(pts)

# check zero mean
sum(∫(ucf)dΩ  )
sum(∫(ucf)dΩg  )

# check compatibility
sum(∫( laplacian(ucf))dΩ  )
sum(∫( laplacian(ucf))dΩg  )

writevtk(Ω,dir*"/poisson",cellfields=["u"=>ucf],append=false)

################################################################################
#### FE spaces on parametric space
################################################################################

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2, constraint=:zeromean)
U = TrialFESpace(V)

T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
S = TrialFESpace(T)

X = MultiFieldFESpace([S,U])
Y = MultiFieldFESpace([T,V])


### Sigma_exact
_sigma_ex = surface_gradient(ucf,m)

writevtk(Ω,dir*"/poisson",cellfields=["ss"=>_sigma_ex],append=false)

### dual form
_rhs = -1.0*surface_divergence(_sigma_ex,m)

biformX((s,u),(t,v)) = (  ∫( s⋅t  )dΩg + ∫( wave_divergence(t,m)*u )dΩ
                            + ∫( surface_divergence(s,m)*v  )dΩg
                      )
liformY((t,v)) = ∫( -(_rhs*v) )dΩg


op = AffineFEOperator(biformX,liformY,X,Y)
sh,uh = solve(LUSolver(),op)

#### Compute errors
e = sum(∫((uh-ucf)⊙(uh-ucf))dΩ)
eg = sum(∫((uh-ucf)⊙(uh-ucf))dΩg)

writevtk(Ω,dir*"/poisson",cellfields=["u"=>ucf,"uh"=>uh,"eu"=>ucf-uh],append=false)
