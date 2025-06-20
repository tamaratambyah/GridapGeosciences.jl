## on a manifold: sphere
function uex(x)
  if x[1] < π
    return x[1]*(π-x[1])
  else
    return (x[1]-π)*(x[1]-2*π)
  end
end
function _uex(x)
  if x[1] < π
    return -x[1]*(π-x[1])
  else
    return -(x[1]-π)*(x[1]-2*π)
  end
end

metric_func(x) = TensorValue{2,2}(1,0.0,0.0,(cos(x[1]))^2 )

model = CartesianDiscreteModel((-π/2,π/2, 0,2*π ), (16,16), isperiodic=(true,true))
Ω = Triangulation(model)
m = Metric(metric_func,Ω)

dΩ = Measure(Ω,degree)
dΩg =  Measure(m,Ω,degree)

uex(x) = cos(x[2])
ucf = CellField(uex,Ω)
# check compatibility
sum( ∫(-ucf)dΩ - ∫( surface_laplacian(ucf,m))dΩ  )

writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex],append=false)


_rhs = -ucf -1.0*(surface_laplacian(ucf,m))

#### FE Problem -- no lagrange multiplers
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(-u*v)dΩg +  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫(  _rhs*v )dΩg
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

A = get_matrix(op)
b = get_vector(op)

evals = eigvals(Array(A))
sum(A*(3*ones(size(b))))
sum(b)

uh = solve(LUSolver(),op)
sum( ∫(uh )dΩ  )

#### Compute errors
e = uh-uex
sum(∫(e⊙e)dΩ)
sum(∫(e⊙e)dΩg)

writevtk(Ω,dir*"/poisson",
        cellfields=["u"=>uex,"uh"=>uh,"e"=>e],append=false)




################################################################################
#### Method 4
#### Mixed form -- with lagrange multiplers
################################################################################
uex(x) = x[1]*(1-x[1])

ucf = CellField(uex,Ω)
# check compatibility
sum( ∫( -ucf)dΩ - ∫( laplacian(ucf))dΩ  )
-1.0*sum( ∫( ucf)dΩ + ∫( laplacian(ucf))dΩ  )
sum( ∫( ucf)dΩ) + sum(∫( laplacian(ucf))dΩ  )

# dual form -- with periodicity, force zero mean
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2)
U = TrialFESpace(V)

T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
S = TrialFESpace(T)

Λ = ConstantFESpace(model)
M = TrialFESpace(Λ)

X = MultiFieldFESpace([S,U,M])
Y = MultiFieldFESpace([T,V,Λ])

### Sigma_exact
sigma_ex(x) = gradient(uex)(x)

_X = MultiFieldFESpace([S,M])
_Y = MultiFieldFESpace([T,Λ])

biformS((s,μ),(t,λ)) = ∫( s⋅t )dΩ + ∫( divergence(s)*λ )dΩ  + ∫( divergence(t)*μ )dΩ
liformS((t,λ)) = ∫( sigma_ex ⋅ t )dΩ
op = AffineFEOperator(biformS,liformS,_X,_Y)
sigma_exh,μh = solve(LUSolver(),op)

# sum(∫((sigma_exh-sigma_ex)⊙(sigma_exh-sigma_ex))dΩ)


### dual form
_rhs = _uex -1.0*divergence(sigma_exh)

biformX((s,u,μ),(t,v,λ)) = ( ∫( u*v )dΩ  + ∫( s⋅t + divergence(t)*u )dΩ
                            + ∫( divergence(s)*v  )dΩ
                            + ∫(v*μ)dΩ + ∫(λ*u)dΩ
                      )
liformY((t,v,λ)) = ∫( -(_rhs*v) )dΩ  + ∫(λ*uex)dΩ

op = AffineFEOperator(biformX,liformY,X,Y)
sh,uh,μh = solve(LUSolver(),op)

e = uh-uex
sum(∫(e⊙e)dΩ)
