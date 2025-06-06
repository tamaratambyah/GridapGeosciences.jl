"""
Consider mapping a square x,y ∈ [0,1]^2 to a plane via a quadratic mapping:
  φ(x,y) = (x,y,x^2 + y^2)
This mapping has a non-constant metric.
Want to test the surface_divergence in the term
   ∫( p*surf_div(v) )dΩg
To avoid unknown sources of error, write the diff operators in the biform out explicitly
Find the evaluation of the biform fails for multifields
"""

using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields, Gridap.TensorValues

model_parametric = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,2)))
Ω_parametric = Triangulation(model_parametric)
dΩ_parametric = Measure(Ω_parametric,2)
pt = Point(0.0,0.0)


### quadratic mapping, φ(x,y) = (x,y,x^2 + y^2)
metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )
sq_meas_func(x) = sqrt(meas(metric_func(x)))

m_cf = CellField(metric_func,Ω_parametric) # 2x2 tensor
sq_meas = CellField(sq_meas_func,Ω_parametric ) # scalar


################################################################################
#### Consider standard cell fields
################################################################################
v(x) = VectorValue(x[1]*x[2],x[2]^2) # vector valued function
_v = map(p->GenericField(v), collect(1:num_cells(model_parametric)))
v_cf = GenericCellField(_v,Ω_parametric,PhysicalDomain())

f = sq_meas*v
f(pt)

divf = divergence(f)
divf(pt)

dc = ∫(divf * sq_meas)dΩ_parametric
sum(dc)

surfdiv = 1/sq_meas * divergence(f)
surfdiv(pt)

dc = ∫(surfdiv * sq_meas)dΩ_parametric
sum(dc)

### seems to be working ...


################################################################################
#### Consider fe spaces
### it fails for multifields, but works for single fields ...
################################################################################
RT = ReferenceFE(raviart_thomas,Float64,1)
DG = ReferenceFE(lagrangian,Float64,1)

V = FESpace(model_parametric, RT, conformity=:Hdiv)
Q = FESpace(model_parametric, DG, conformity=:L2)

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


biform1((u,p),(v,q)) = ∫( (p*( (1/sq_meas * divergence(sq_meas * v) ) )) *sq_meas )dΩ_parametric
biform2((u,p),(v,q)) = ∫( (p*( (1/sq_meas * divergence( v) ) )) *sq_meas )dΩ_parametric ## removed product in divergence

### single fields
dc1 = ∫(  p_trial * (1/sq_meas * divergence(sq_meas * v_test) )* sq_meas )dΩ_parametric
biform1((u_trial,p_trial),(v_test,q_test))# works
biform2((u_trial,p_trial),(v_test,q_test)) # works

### multifield
biform1(x_trial,y_test) # fails
biform2(x_trial,y_test) # works -- removed product in divergence
