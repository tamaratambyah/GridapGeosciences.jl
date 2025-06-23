using Gridap
using Plots, LaTeXStrings
using DrWatson
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns


n = 8

# Domain initialisation -- 2D
domain = (0,1,0,1)
partition = (n,n)

uex(x) = VectorValue(0.0,x[2]*(1-x[2]))
pex(x) = 1.0 +  0.01*x[1]*(1-x[1])


##
model = CartesianDiscreteModel(domain, partition)

Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ,degree)
n_Γ = get_normal_vector(Γ)

ucf = CellField(uex,Ω)
pcf = CellField(pex,Ω)

writevtk(Ω,dir*"/wave",cellfields=["u"=>uex,"p"=>pex],append=false)


u0 = ucf + gradient(pcf)
p0 = pcf + divergence(ucf)

### FE problem
V = TestFESpace(model,ReferenceFE(raviart_thomas,Float64,p),conformity=:HDiv)
U = TrialFESpace(V)

Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,p),conformity=:L2)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V,Q])

# Weak formulation for u
wave_biform((u,p),(v,q)) = ( ∫( u⋅v )dΩ
                   + ∫( -1.0*p*(divergence(v))   )dΩ
                  + ∫( p*q )dΩ
                  + ∫( divergence(u)*q )dΩ
                  )
wave_liform((v,q)) = ∫( u0⋅v + p0*q  )dΩ -  ∫( pcf*(v⋅n_Γ) )dΓ

op = AffineFEOperator(wave_biform,wave_liform,X,Y)
uh, ph = solve(LUSolver(),op)

# Error
eu =  sum(∫((uh-uex)⊙(uh-uex))dΩ)
ep = sum(∫((ph-pex)⊙(ph-pex))dΩ)


writevtk(Ω,dir*"/wave",
    cellfields=["u"=>uex,"p"=>pex,
                "uh"=>uh, "ph"=>ph,
                "eu"=>uh-uex,"ep"=>ph-pex],append=false)
