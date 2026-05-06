"""
Test the overloads that are triggered by a specific combination of the product
rule in the divergence
"""

module OverloadTests


using Gridap
using Test

p_fe = 1
model = CartesianDiscreteModel((0,1,0,1),(2,2))
Ω = Triangulation(model)
dΩ = Measure(Ω,2*p_fe)

test_func(x) = x[1]
cf = CellField(test_func,Ω)

### Finite element spaces
RT = ReferenceFE(raviart_thomas,Float64,1)
DG = ReferenceFE(lagrangian,Float64,1)

V = FESpace(model, RT, conformity=:Hdiv)
Q = FESpace(model, DG, conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U, P])
Y = MultiFieldFESpace([V, Q])

u_trial = get_trial_fe_basis(U)
v_test = get_fe_basis(V)
p_trial = get_trial_fe_basis(P)
q_test = get_fe_basis(Q)

x_trial = get_trial_fe_basis(X)
y_test = get_fe_basis(Y)

### Consider biforms that contain the product rule in the divergence
biform1((u,p),(v,q)) = ∫((p*( ( divergence(cf * v) ) )) )dΩ
biform2((u,p),(v,q)) = ∫((q*( ( divergence(cf * u) ) )) )dΩ # Transpose

### Evaluate the biforms with singlefields
biform1((u_trial,p_trial),(v_test,q_test))# works
biform2((u_trial,p_trial),(v_test,q_test))# works

### Evaluate the birforms with mulitfields -- these are broken
@test_skip( @test_broken biform1(x_trial,y_test) )
@test_skip( @test_broken biform2(x_trial,y_test) )

### Now include the overloads in GridapGeosciences, to evluate with multifield
using GridapGeosciences
biform1(x_trial,y_test)
biform2(x_trial,y_test)

@test true



end # module
