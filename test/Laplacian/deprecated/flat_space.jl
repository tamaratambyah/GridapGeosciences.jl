using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using Gridap.ReferenceFEs, Gridap.Polynomials, Gridap.CellData

using GridapGeosciences


include("../convergence_tools.jl")

dir = datadir("HodgeLaplacian_flat")
!isdir(dir) && mkdir(dir)

function dx(model::UnstructuredDiscreteModel{3,3})
  1/( num_cells(model)^(1/3) )
end

model0 = UnstructuredDiscreteModel( CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (2, 2, 2)))
model1 = UnstructuredDiscreteModel( CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (4, 4, 4)))
model2 = UnstructuredDiscreteModel( CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (16, 16, 16)))

models = [model2, model1, model0]
ps = [0,1]
ls = LUSolver()

function uX(x)
  VectorValue(x[1]^3,x[2]^2*x[1],x[1]*x[3]^2)
  # VectorValue(x[1]^3,cos(2*π*x[1])*sin(2*π*x[1]),sin(2*π*x[1]))
end


p_convergence_test([true],ps,models,flat_hodge_laplacian,dir,uX,ls,true)
plot_convergence_from_saved(dir,"convergence",["u:", "s:"])

function flat_hodge_laplacian(
  model,
  p_fe::Int,dir::String,uX::Function,ls=LUSolver(),return_vtk=false)

  degree = 4*(p_fe + 1)
  if p_fe == 0
    degree = 8
  end
  @check degree > 0 "Zero quad!!"

  ## finite element solver
  Ω = Triangulation(model)
  dΩ = Measure(Ω,2*degree)
  Ω_error = Triangulation(model)
  dΩ_error = Measure(Ω_error,4*degree)

  tags = ["boundary"]
  Γ = BoundaryTriangulation(model,tags=tags)
  dΓ = Measure(Γ,2*degree)
  nΓ = get_normal_vector(Γ)


  divuX(x) = (∇⋅uX)(x)
  curluX(x) = (∇×uX)(x)
  ccurl(x) = (∇×(∇×uX))(x)
  rhs(x) = ccurl(x) - gradient(divuX)(x)

  ## cellfields
  rhs_cf = CellField(rhs,Ω)
  u_cf = CellField(uX,Ω)
  div_cf =  CellField(divuX,Ω)
  sigma_cf = -div_cf
  curlu_cf = CellField(curluX,Ω)


  ## FE spaces
  T = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  S = TrialFESpace(T)

  R = TestFESpace(Ω, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl)
  H = TrialFESpace(R)


  ### Multifield
  X = MultiFieldFESpace([S,H])
  Y = MultiFieldFESpace([T,R])

  biform_x((s,u),(t,v)) = (
                  ∫( (s*t)  )dΩ
                - ∫( ∇(t)⋅u  )dΩ
                + ∫( curl(u)⋅curl(v) )dΩ
                + ∫( gradient(s)⋅v )dΩ
                  )
  liform_x((t,v)) = (
                ∫( rhs_cf⋅v  )dΩ
                + ∫( v⋅( curlu_cf ×nΓ) )dΓ
                - ∫(( t*(u_cf⋅nΓ) )  )dΓ
                  )


  op = AffineFEOperator(biform_x,liform_x,X,Y)
  sh,uh = solve(ls,op)


  _e = sigma_cf - sh
  el2_s = sqrt(sum(∫( (_e*_e)  )dΩ_error))

  _e = uh - u_cf
  el2_u = sqrt(sum(∫( (_e⋅_e)  )dΩ_error))

  if return_vtk
    cellfields =  ["u"=>u_cf,
    "uh"=>uh,
    "eu"=>uh-u_cf,
    "sh"=>sh, "s"=>sigma_cf, "e"=>sh-sigma_cf,
                  ]
    writevtk(Ω,dir*"/flat_model_nref",
            cellfields=cellfields,
            append=false)
  end

  return el2_u, el2_s, false

end
