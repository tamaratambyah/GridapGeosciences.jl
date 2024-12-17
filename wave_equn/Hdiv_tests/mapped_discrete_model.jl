using Gridap
using GridapSolvers
using GridapGeosciences
using LinearAlgebra
using DrWatson


out_dir = datadir("Hdiv")

nc = 4
order = 0
degree = 6

domain = (-1,1,-1,1,-1,1)
cells  = (nc,nc,nc)
cube  = CartesianDiscreteModel(domain,cells)



u_exact = VectorValue(1,2,3)

function f(xyz) # map to sphere
  x,y,z = xyz
  xₛ = x*sqrt(1.0-y^2/2-z^2/2+y^2*z^2/(3.0))
  yₛ = y*sqrt(1.0-z^2/2-x^2/2+x^2*z^2/(3.0))
  zₛ = z*sqrt(1.0-x^2/2-y^2/2+x^2*y^2/(3.0))
  Point(xₛ,yₛ,zₛ)
end


new_model = Gridap.Geometry.MappedDiscreteModel(cube,f)


Ω = Triangulation(new_model)
dΩ = Measure(Ω, degree)
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

u_cf = CellField(u_exact,Ω)


rt_reffe = ReferenceFE(raviart_thomas, Float64, 0)

V = FESpace(new_model, rt_reffe; conformity=:Hdiv)
U = TrialFESpace(V)


# interpolate
u_int = interpolate(u_exact,U)

# project onto Hdiv
a(u,v) = ∫( u⋅v )dΩ
function l(f,dΩ)
 (v) -> ∫( f⋅v  )dΩ
end

op = AffineFEOperator(a,l(u_exact,dΩ),U,V)
u_proj = solve(LUSolver(),op)


writevtk(Ω,joinpath(out_dir,"mapped_sphere"),
                cellfields=["u_cf"=>u_cf,"u_int"=>u_int,"u_proj"=>u_proj,
                            "err_proj"=> u_cf-u_proj,
                            "err_int"=> u_cf-u_int,
                           ],append=false)
