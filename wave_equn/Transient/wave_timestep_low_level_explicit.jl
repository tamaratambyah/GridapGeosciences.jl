"""
low level time stepping using SSPRK2 - explicit methods
	U_{n+1} = (I + Δt B) Uₙ

  ϕ₁ = ϕₙ + Δt B(ϕₙ)
  ϕ₂ = ϕ₁ + Δt B(ϕ₁)
  ϕ_{n+1} = 0.5 ( ϕₙ + ϕ₂ )

"""

using Gridap
using GridapSolvers
using GridapGeosciences
using LinearAlgebra
using DrWatson

include("../createpvd.jl")



p0(t) = x -> 1.0 + 0.0001*exp(-5*((1-x[1])^2+(-x[2])^2+(-x[3])^2))
u0(t) = x -> VectorValue(0.0 , 0.0, 0.0 )

n = 16
p = 1
degree = 6# 2*(p+1)

Nstep = 2000
dt = π/Nstep

out_dir = datadir("wave_transient_n$(n)_p$(p)_vanka")
!isdir(out_dir) && mkdir(out_dir)
pvd = createpvd(out_dir)

model = CubedSphereDiscreteModel(n)
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)
dω = Measure(Ω,degree,ReferenceDomain())
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

rt_reffe = ReferenceFE(raviart_thomas, Float64, p)
lg_reffe = ReferenceFE(lagrangian, Float64, p)

V = FESpace(model, rt_reffe, conformity=:Hdiv)
U = TrialFESpace(V)

Q = FESpace(model, lg_reffe; conformity=:L2)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U, P])
Y = MultiFieldFESpace([V, Q])

biform((u,p),(v,q)) = ∫( u⋅v + p*q  )dΩ
liform((v,q),(un,pn)) = ∫( un⋅v + pn*q  )dΩ + dt*( ∫( (DIV(v))*pn - (DIV(un))*q )dω )

### initial conditions
xh0 = interpolate_everywhere([u0(0.0),p0(0.0)], X)
uh0,ph0 = xh0

### time stepper set up
xn = get_free_dof_values(xh0)

x = get_trial_fe_basis(X)
y = get_fe_basis(Y)
assem = SparseMatrixAssembler(X,Y)
matdata = Gridap.FESpaces.collect_cell_matrix(X, Y, biform(x, y))
vecdata = Gridap.FESpaces.collect_cell_vector(Y,liform(y,xh0))
A = Gridap.FESpaces.assemble_matrix(assem, matdata)
b = Gridap.FESpaces.allocate_vector(assem,vecdata)

# work vectors
x1 = Gridap.FESpaces.allocate_in_domain(A); fill!(x1,0.0)
x2 = Gridap.FESpaces.allocate_in_domain(A); fill!(x2,0.0)

## solver
# PD = PatchDecomposition(model)
# P = GridapSolvers.VankaSolver(X,PD)

P = JacobiLinearSolver() #GridapSolvers.VankaSolver(Y)
ls = GMRESSolver(20;Pl=P,maxiter=1000,atol=1e-14,rtol=1.e-14,verbose=1)
# ls = LUSolver()


# facotrise matrix
ss = symbolic_setup(ls,A)
ns = numerical_setup(ss,A)

@time for t = 1:Nstep
  println(t)

  xhn = FEFunction(X,xn)
  uh,ph = xhn
  pvd[t-1] = createvtk(Ω,joinpath(out_dir,"wave_t$(t-1)"),
                  cellfields=["uh"=>uh,"ph"=>ph],append=false)


  ### stage 1:
  # Gridap.FESpaces.assemble_matrix!(biform(x,y), A,assem,X,Y)
  Gridap.FESpaces.assemble_vector!(liform(y,xhn), b,assem,Y)
  solve!(x1,ns,b)
  x1h = FEFunction(X,x1)

  ### stage 2:
  # Gridap.FESpaces.assemble_matrix!(biform(x,y), A,assem,X,Y)
  Gridap.FESpaces.assemble_vector!(liform(y,x1h), b,assem,Y)
  solve!(x2,ns,b)

  ### combine in place stage solutions
  xn .= 0.5 .* ( xn .+ x2 )

end

xhn = FEFunction(X,xn)
uh,ph = xhn
pvd[Nstep] = createvtk(Ω,joinpath(out_dir,"wave_t$(Nstep)"),
                cellfields=["uh"=>uh,"ph"=>ph],append=false)


make_pvd(out_dir,"wave_t","output",1)
