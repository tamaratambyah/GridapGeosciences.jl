using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity
using LinearAlgebra
global RADIUS = 1.0

function cube_to_alphabeta(XYZ)
  α,β = -10.0,-10.0
  if p == 1
  elseif p == 3
  end
  return Point(α,β)
end

function forward_map(αβ,p)
  α,β = αβ
  X,Y,Z = 0.0,0.0,0.0

  rho = sqrt( 1 + (tan(α))^2 + (tan(β))^2  )

  if p == 1
    X = 1/rho
    Y = X * tan(α)
    Z = X * tan(β)
  elseif p == 3
    Y = 1/rho
    X = -Y * tan(α)
    Z = Y * tan(β)
  elseif p == 4
    X = -1/rho
    Y = -X * tan(β)
    Z = -X * tan(α)
  elseif p == 6
    Y = -1/rho
    X = Y * tan(β)
    Z = -Y * tan(α)
  end
  XYZ = RADIUS*Point(X,Y,Z)
  return XYZ

end

function forward_jacobian(αβ,p)
  α,β = αβ
  dXda,dYda,dZda = 0.0,0.0,0.0
  dXdb,dYdb,dZdb = 0.0,0.0,0.0

  rho = sqrt( 1 + (tan(α))^2 + (tan(β))^2  )

  dda = -1.0*tan(α)*(sec(α))^2 / rho^3
  ddb = -1.0*tan(β)*(sec(β))^2 / rho^3
  factor = 1/rho

  if p == 1
    dXda = dda
    dXdb = ddb
    dYda = dda*tan(α) + factor*(sec(α))^2
    dYdb = ddb*tan(α)
    dZda = dda*tan(β)
    dZdb = ddb*tan(β) + factor*(sec(β))^2
  elseif p == 3
    dYda = dda
    dYdb = ddb
    dXda = -1.0*( dda*tan(α) + factor*(sec(α))^2 )
    dXdb = -1.0*( ddb*tan(α) )
    dZda = dda*tan(β)
    dZdb = ddb*tan(β) + factor*(sec(β))^2
  elseif p == 4
    dXda = -dda
    dXdb = -ddb
    dYda = -( -dda*tan(β) )
    dYdb = -( -ddb*tan(β) + -factor*(sec(β))^2 )
    dZda = -( -dda*tan(α) + -factor*(sec(α))^2 )
    dZdb = -( -ddb*tan(α) )
  elseif p == 6
    dYda = -dda
    dYdb = -ddb
    dXda = -dda*tan(β)
    dXdb = -ddb*tan(β) + -factor*(sec(β))^2
    dZda = -( -dda*tan(α) + -factor*(sec(α))^2 )
    dZdb = -( -ddb*tan(α) )
  end

  cov_basis_1 = RADIUS*Point(dXda,dYda,dZda)
  cov_basis_2 = RADIUS*Point(dXdb,dYdb,dZdb)
  return cov_basis_1,cov_basis_2

end


function metric(αβ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
  E = cov_basis_1 ⋅ cov_basis_1
  F = cov_basis_1 ⋅ cov_basis_2
  G = cov_basis_2 ⋅ cov_basis_2

  g = TensorValue{2,2,Float64}(E,F,F,G)

  return g
end

function meas_metric(αβ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
  E = cov_basis_1 ⋅ cov_basis_1
  F = cov_basis_1 ⋅ cov_basis_2
  G = cov_basis_2 ⋅ cov_basis_2

  sqrt_detg = sqrt(E*G - F*F)

  return sqrt_detg
end

function inv_metric(αβ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
  E = cov_basis_1 ⋅ cov_basis_1
  F = cov_basis_1 ⋅ cov_basis_2
  G = cov_basis_2 ⋅ cov_basis_2

  detg = (E*G - F*F)

  ginv = 1/detg * TensorValue{2,2}(G,-F,-F,E)

  return ginv
end





function inverse_map(XYZ,p)
  X,Y,Z = XYZ
  α,β = -10.0,-10.0

  if p == 1
    α = atan(Y,X)
    β = atan(Z,X)
  elseif p == 3
    α = atan(-X,Y)
    β = atan(Z,Y)
  elseif p == 4
    α = atan(-Z,X)
    β = atan(-Y,X)
  elseif p == 6
    α = atan(-Z,Y)
    β = atan(X,Y)
  end

  αβ = Point(α,β)
  return αβ


end

function inverse_jacobian(XYZ,p)
  X,Y,Z = XYZ
  dadX,dadY,dadZ = 0.0,0.0,0.0
  dbdX,dbdY,dbdZ = 0.0,0.0,0.0

  if p == 1
    dadX = -1.0* Y / (X^2 + Y^2)
    dadY = X / (X^2 + Y^2)
    dadZ = 0.0
    dbdX = -1.0 * Z / (X^2 + Z^2)
    dbdY = 0.0
    dbdZ = X / (Z^2 + X^2)
  elseif p == 3
    dadX = -1.0* Y / (X^2 + Y^2)
    dadY = X / (X^2 + Y^2)
    dadZ = 0.0
    dbdX = 0.0
    dbdY = -1.0*Z/( Z^2 + Y^2)
    dbdZ = Y / (Z^2 + Y^2)
  elseif p == 4
    dadX = Z / (X^2 + Z^2)
    dadY = 0.0
    dadZ = -1.0*X / (X^2 + Z^2)
    dbdX = Y / (X^2 + Y^2)
    dbdY = -1.0*X/( X^2 + Y^2)
    dbdZ = 0.0
  elseif p == 6
    dadX = 0.0
    dadY = Z / (Y^2 + Z^2)
    dadZ = -1.0*Y / (Y^2 + Z^2)
    dbdX = Y / (X^2 + Y^2)
    dbdY = -1.0*X/( X^2 + Y^2)
    dbdZ = 0.0
  end

  contra_basis_1 = Point(dadX,dadY,dadZ)
  contra_basis_2 = Point(dbdX,dbdY,dbdZ)
  return contra_basis_1,contra_basis_2
end

function convarient_components(vec_phys,αβ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)
  u_cov = vec_phys ⋅ cov_basis_1
  v_cov = vec_phys ⋅ cov_basis_2
  return u_cov,v_cov
end

function contrvarient_components(vec_phys,XYZ,p)
  contra_basis_1,contra_basis_2 = inverse_jacobian(XYZ,p)
  u_contra = vec_phys ⋅ contra_basis_1
  v_contra = vec_phys ⋅ contra_basis_2
  return u_contra,v_contra
end


function cov_to_phys(vec_phys,αβ,p)
  XYZ = forward_map(αβ,p)
  u_cov,v_cov = convarient_components(vec_phys,αβ,p)

  contra_basis_1,contra_basis_2 = inverse_jacobian(XYZ,p)
  u = u_cov*contra_basis_1[1] + v_cov*contra_basis_2[1]
  v = u_cov*contra_basis_1[2] + v_cov*contra_basis_2[2]
  w = u_cov*contra_basis_1[3] + v_cov*contra_basis_2[3]

  return Point(u,v,w)

end

function contra_to_phys(vec_phys,αβ,p)
  XYZ = forward_map(αβ,p)
  u_contra,v_contra = contrvarient_components(vec_phys,XYZ,p)
  cov_basis_1,cov_basis_2 = forward_jacobian(αβ,p)

  u = u_contra*cov_basis_1[1] + v_contra*cov_basis_2[1]
  v = u_contra*cov_basis_1[2] + v_contra*cov_basis_2[2]
  w = u_contra*cov_basis_1[3] + v_contra*cov_basis_2[3]

  return Point(u,v,w)
end



p = 1
αβ = Point(-π/4,π/4)

XYZ = forward_map(αβ,p)
cov_v1,cov_v2 = forward_jacobian(αβ,p)

αβ = inverse_map(XYZ,p)
con_v1,con_v2 = inverse_jacobian(XYZ,p)


######## mapping of tangent vector. At node 4
vec = Point(-2.0,1.0,1.0)

# panel 1:
p = 1
αβ = Point(π/4,π/4)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 3:
p = 3
αβ = Point(-π/4,π/4)
cov_to_phys(vec,αβ,p)

contra_to_phys(vec,αβ,p)

######## mapping of tangent vector. At node 3
vec = Point(1.0,2.0,1.0)

# panel 1:
p = 1
αβ = Point(-π/4,π/4)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 6:
p = 6
αβ = Point(π/4,-π/4)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)


######## mapping of tangent vector. At node 5
vec = Point(1.0,1.0,2.0)

# panel 4:
p = 4
αβ = Point(π/4,-π/4)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

# panel 6:
p = 6
αβ = Point(π/4,π/4)
cov_to_phys(vec,αβ,p)
contra_to_phys(vec,αβ,p)

###############################################################################
u(X) = X[1]*X[2]*X[3]

################################################################################
##### Single panel: 1, 3
################################################################################
## force the panel number below
function uex(αβ)
  XYZ = forward_map(αβ,3)
  u(XYZ)
end

function panel1_meas_metric(αβ)
  meas_metric(αβ,3)
end

function panel1_inv_metric(αβ)
  inv_metric(αβ,3)
end

single_panel_model = CartesianDiscreteModel(( -π/4,π/4,-π/4,π/4),(20,20))

p = 2
degree = 4*p

Ω = Triangulation(single_panel_model)
dΩ = Measure(Ω,degree)
pts = get_cell_points(Ω)

ucf = CellField(uex,Ω)
meas_metric_cf = CellField(panel1_meas_metric,Ω)
inv_metric_cf = CellField(panel1_inv_metric,Ω)

surflap = 1/meas_metric_cf * divergence( meas_metric_cf *( inv_metric_cf ⋅ gradient(ucf) ) )

sum(∫( ucf*meas_metric_cf + surflap*meas_metric_cf  )dΩ)

rhs = ucf + surflap


V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(u*v*meas_metric_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_metric_cf )dΩ
poisson_liform(v) = ∫(  (rhs*v)*meas_metric_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

uh = solve(LUSolver(),op)


l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-ucf,dΩ)



################################################################################
##### Single panel: 4, 6
################################################################################
single_panel_model = UnstructuredDiscreteModel( CartesianDiscreteModel(( -π/4,π/4,-π/4,π/4),(20,20)) )
grid = get_grid(single_panel_model)
cmaps = get_cell_map(grid)
nodes = Geometry.get_node_coordinates(grid)

A = [0 1
    1 0]

swap = fill(SwapField(TensorValue(A)),length(cmaps))
new_cmaps = lazy_map(∘,swap,cmaps)

evaluate(cmaps[2],Point(0,0))
evaluate(new_cmaps[2],Point(0,0))

new_cell_coords = lazy_map(evaluate,new_cmaps,get_cell_ref_coordinates(grid))./1
new_nodes = map(x->Point(x[2],x[1]),nodes)
new_grid = Geometry.UnstructuredGrid(new_nodes,get_cell_node_ids(grid),get_reffes(grid),get_cell_type(grid),OrientationStyle(grid),
                    nothing,new_cmaps)


new_single_panel_model = UnstructuredDiscreteModel(new_grid,get_grid_topology(single_panel_model),get_face_labeling(single_panel_model))

## force the panel number below
function uex(αβ)
  XYZ = forward_map(αβ,6)
  u(XYZ)
end

function panel1_meas_metric(αβ)
  meas_metric(αβ,6)
end

function panel1_inv_metric(αβ)
  inv_metric(αβ,6)
end


p = 2
degree = 4*p

Ω = Triangulation(new_single_panel_model)
dΩ = Measure(Ω,degree)
pts = get_cell_points(Ω)

ucf = CellField(uex,Ω)
meas_metric_cf = CellField(panel1_meas_metric,Ω)
inv_metric_cf = CellField(panel1_inv_metric,Ω)

surflap = 1/meas_metric_cf * divergence( meas_metric_cf *( inv_metric_cf ⋅(gradient(ucf)) ) )

sum(∫( ucf*meas_metric_cf + surflap*meas_metric_cf  )dΩ)

rhs = ucf + surflap

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(u*v*meas_metric_cf)dΩ -  ∫( ( ∇(v)⋅ (inv_metric_cf⋅ ∇(u))  )*meas_metric_cf )dΩ
poisson_liform(v) = ∫(  (rhs*v)*meas_metric_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

uh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-ucf,dΩ)
