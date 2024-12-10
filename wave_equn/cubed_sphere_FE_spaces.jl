"""
Compares CellFields evaluations and interpolation of functions on cubed sphere
mesh
"""

using Gridap
using GridapGeosciences
using DrWatson

function p_exact_solution(p)
  if p == 0
    return x -> 1.0
  # elseif p == 1
  #   return x -> 1.0 + 0.1*x[1]
  elseif p == 2
    return  x -> 1.0 + 0.1*x[1]*(1-x[1])
  else
    println("non periodic p")
  end
end


function u_exact_solution(p,d,model)
  if d == 2
    @assert typeof(model) <: CartesianDiscreteModel
    _u_exact_2D(p)
  else
    @assert d == 3
    _u_exact_3D(p)
  end
end


function _u_exact_2D(p)
  if p == 0
    return x -> VectorValue(1.0 , 2.0 )
  # elseif p == 1
  #   return x -> VectorValue(x[1] , 0.0 )
  elseif p == 2
    return x -> VectorValue(x[1]*(1-x[1]), 0.0 )
  else
    println("non periodic p")
  end
end

function _u_exact_3D(p)
  if p == 0
    return x -> VectorValue(1.0 , 2.0, 3.0 )
  # elseif p == 1
  #   return x -> VectorValue(x[1] , 0.0, 0.0 )
  elseif p == 2
    # return x -> VectorValue(x[1]*(1-x[1]) , 0.0, 0.0 )
    return x -> VectorValue(x[1]*x[1] , 0.0, 0.0 )
    # return x -> VectorValue(x[2]x[3]- 3x[1]x[1]x[2]x[3] / (x[1]x[1] + x[2]x[2] + x[3]x[3]),
    #                         x[1]x[3] - 3x[1]x[2]x[2]x[3]/ (x[1]x[1] + x[2]x[2] + x[3]x[3]),
    #                         x[1]x[2] - 3x[1]x[2]x[3]x[3] / (x[1]x[1] + x[2]x[2] + x[3]x[3]) )

  else
    println("non periodic p")
  end
end


function get_cart_model(n,d)
  if d == 2
    return CartesianDiscreteModel((0,1,0,1), (n,n),isperiodic=(true,true))
  else
    @assert d == 3
    return  CartesianDiscreteModel((0,1,0,1,0,1), (n,n,n),isperiodic=(true,true,true))
  end
end

function generate_output(cfs,model)

  if typeof(model) <: GridapGeosciences.AnalyticalMapCubedSphereDiscreteModel
    name = "cubed_sphere_wave_exact"
  else
    name = "cart_wave_exact"
  end

  Ω = Triangulation(model)
  writevtk(Ω,joinpath(datadir("models"),name),
    cellfields=cfs,
                append=false)
end

l2(e,dΩ) = sum(∫(e⊙e)dΩ)

n = 4
p = 2
d = 3

# model = get_cart_model(n,d)
model = CubedSphereDiscreteModel(n)

p_exact = p_exact_solution(p)
u_exact = u_exact_solution(p,d,model)

degree = 2*(p+1)

Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

rt_reffe = ReferenceFE(raviart_thomas, Float64, p)
lg_reffe = ReferenceFE(lagrangian, Float64, p)

V = FESpace(model, rt_reffe; conformity=:Hdiv)
U = TrialFESpace(V)

Q = FESpace(model, lg_reffe; conformity=:L2)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V,Q])

# project solution into Hdiv
a0((u,p),(v,q)) = ∫( u⋅v + p*q )dΩ
l0((v,q)) = ∫( u_exact⋅v + p_exact*q )dΩ
op = AffineFEOperator(a0,l0,X,Y)
uh, ph = solve(LUSolver(),op)

cfs=["u"=>CellField(u_exact,Ω),"u_ex"=>interpolate(u_exact,U),"uh"=>uh,
    "eu"=>CellField(u_exact,Ω)-interpolate(u_exact,U),
    "p"=>CellField(p_exact,Ω),"p_ex"=>interpolate(p_exact,P),"ph"=>ph,
    "ep"=>CellField(p_exact,Ω)-interpolate(p_exact,P)]

generate_output(cfs,model)
l2(CellField(u_exact,Ω)-interpolate(u_exact,U),dΩ)
