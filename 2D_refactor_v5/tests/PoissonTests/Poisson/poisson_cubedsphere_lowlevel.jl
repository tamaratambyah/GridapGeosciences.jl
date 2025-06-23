using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../../../src/initialise.jl")
include("../poisson_helpers.jl")



#### Analytic solution with zero mean
p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

n = ns[3]

global RADIUS = 1.0*sqrt(3.0)


################################################################################
#### Analytic parametric space for a single panel
################################################################################
n = 16
model = CartesianDiscreteModel((-π/4,π/4, -π/4,π/4), (n,n))
Ω = Triangulation(model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ,degree)
dΓg = Measure(m,Γ,degree)
n_Γ = get_normal_vector(Γ)

function uex(x)
  if x[1] < 0.0
    return -x[1]*(x[1] + π/4)
  else
    return x[1]*(x[1] - π/4)
  end

end


writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex],append=false)

ucf = CellField(uex,Ω)

# check zero mean and compatibility
sum(∫(ucf)dΩ  ), sum(∫( surface_laplacian(ucf,m))dΩ  )
sum(∫(ucf)dΩg  ), sum(∫( surface_laplacian(ucf,m))dΩg  )


_rhs = -1.0*surface_laplacian(ucf,m)
h = surface_gradient(ucf,m)⋅n_Γ


## zero mean in FE space
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1, constraint=:zeromean)
U = TrialFESpace(V)

poisson_biform(u,v) =  ∫( surface_gradient(u,m)⋅gradient(v) )dΩg
poisson_liform(v) =  ∫( (_rhs*v) )dΩg + ∫( v*h )dΓg

op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
# A = get_matrix(op)
# b = get_vector(op)
# eigvals(Array(A))
# _vec = A*ones(size(b))
# norm(A*ones(size(b)))
# sum(b)

uh = solve(LUSolver(),op)


e = l2(uh-ucf,dΩ)
eg = l2(uh-ucf,dΩg)
writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex,"uh"=>uh,"e"=>uex-uh],append=false)
