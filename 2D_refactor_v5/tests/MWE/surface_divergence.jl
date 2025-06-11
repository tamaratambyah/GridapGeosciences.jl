"""
Consider mapping a square x,y ∈ [0,1]^2 to a plane via a quadratic mapping:
  φ(x,y) = (x,y,x^2 + y^2)
This mapping has a non-constant metric.
Want to test the surface_divergence in the term
   ∫( p*surf_div(v) )dΩg
To avoid unknown sources of error, write the diff operators in the biform out explicitly
Find the evaluation of the biform no longer fails for multifields
"""

using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields, Gridap.TensorValues

################################################################################
#### I think the following extra overloads should eventually go to Gridap's core
#### They are needed to have the product rule working whenever there is the
#### specific combination of Field and VectorBlock{<:Field} that is triggered below
################################################################################
function Gridap.Fields.return_value(k::Broadcasting{<:Operation},
                                    f1::Field,f2::VectorBlock{A},g1::Field,g2::VectorBlock{B}) where {A,B}
  @assert length(f2.array) == length(g2.array)                       
  @assert f2.touched == g2.touched
  
  f2i = Gridap.Fields.testitem(f2)
  g2i = Gridap.Fields.testitem(g2)
  f1f2ig1g2i = Gridap.Fields.return_value(k,f1,f2i,g1,g2i)
  o = Vector{typeof(f1f2ig1g2i)}(undef,size(f2.array))
  for i in eachindex(f2.array)
    if f2.touched[i]
      o[i] = Gridap.Fields.return_value(k,f1,f2.array[i],g1,g2.array[i])
    end
  end
  ArrayBlock(o,f2.touched)
end

function Gridap.Fields.return_cache(k::Broadcasting{<:Operation},
                                    f1::Field,f2::VectorBlock{A},g1::Field,g2::VectorBlock{B}) where {A,B}
  @assert length(f2.array) == length(g2.array)                       
  @assert f2.touched == g2.touched
  
  f2i = Gridap.Fields.testitem(f2)
  g2i = Gridap.Fields.testitem(g2)

  cf1f2ig1g2i = Gridap.Fields.return_cache(k,f1,f2i,g1,g2i)
  f1f2ig1g2i = Gridap.Fields.evaluate!(cf1f2ig1g2i,k, f1,f2i,g1,g2i)
  
  l = Vector{typeof(cf1f2ig1g2i)}(undef,size(f2.array))
  o = Vector{typeof(f1f2ig1g2i)}(undef,size(f2.array))
  for i in eachindex(f2.array)
    if f2.touched[i]
      l[i] = return_cache(k,f1,f2.array[i],g2,g2.array[i])
    end
  end
  ArrayBlock(o,f2.touched),l
end

function Gridap.Fields.evaluate!(cache,k::Broadcasting{<:Operation},
                                 f1::Field,f2::VectorBlock{A},g1::Field,g2::VectorBlock{B}) where {A,B}
  
  o,l = cache
  @assert o.touched == f2.touched
  for i in eachindex(f2.array)
     if f2.touched[i]
       o.array[i] = evaluate!(l[i],k,f1, f2.array[i],g1,g2.array[i])
     end
  end
  o
end
################################################################################
#### End of extra overloads
################################################################################


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


biform1((u,p),(v,q)) = ∫((p*( (1/sq_meas * divergence(sq_meas * v) ) )) *sq_meas)dΩ_parametric
biform2((u,p),(v,q)) = ∫((p*( (1/sq_meas * divergence( v) ) )) *sq_meas)dΩ_parametric  ## removed product in divergence

### single fields
dc1 = ∫(  p_trial * (1/sq_meas * divergence(sq_meas * v_test) )* sq_meas )dΩ_parametric
biform1_scalar(p,v) = ∫(p * (1/sq_meas * divergence(sq_meas * v) )* sq_meas )dΩ_parametric
biform1((u_trial,p_trial),(v_test,q_test))# works
biform2((u_trial,p_trial),(v_test,q_test)) # works

### multifield
A1=assemble_matrix(biform1,X,Y) # already works
A2=assemble_matrix(biform1_scalar,P,V)
@assert A1[1:40,41:end]≈A2

biform2(x_trial,y_test) # works -- removed product in divergence


# DEBUG statements (uncomment to use)
# v_test_y_test, _ = y_test  
# x_x = map(p->GenericField(x->x[1]), collect(1:num_cells(model_parametric)))
# gcf = GenericCellField(x_x, Ω_parametric, ReferenceDomain())
# gradient(gcf)
# gradient(1.0*v_test_y_test)
# prod = 2.0*v_test_y_test
# prod_data = Gridap.CellData.get_data(prod)
# lazy_map(Broadcasting(gradient),prod_data)
# a=prod_data
# f = a.args
# g = map(i->lazy_map(Broadcasting(∇),i),f)
# k(F1,F2,G1,G2) = Gridap.Fields.product_rule(*,F1,F2,G1,G2)
# fi = map(Gridap.Fields.testitem,(f[1],f[2],g[1],g[2]))
# T = return_cache(Broadcasting(Operation(k)), (fi[1],fi[2].array[1],fi[3],fi[4].array[1])...)
# evaluate!(T, Broadcasting(Operation(k)), fi[1],fi[2].array[1],fi[3],fi[4].array[1])