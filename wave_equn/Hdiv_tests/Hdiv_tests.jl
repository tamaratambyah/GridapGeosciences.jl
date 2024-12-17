using Gridap
using GridapSolvers
using GridapGeosciences
using LinearAlgebra
using DrWatson

# include("TensorFESpaces.jl")
# include("TensorQuadratures.jl")

out_dir = datadir("Hdiv")

n = 4
order = 0
degree = 6

# model = CartesianDiscreteModel((0,1,0,1,0,1), (n,n,n))
model = CubedSphereDiscreteModel(n)
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

rt_reffe = ReferenceFE(raviart_thomas, Float64, order)

V = FESpace(model, rt_reffe; conformity=:Hdiv)
U = TrialFESpace(V)

u_exact(x) = VectorValue(1.0, 2.0, 3.0)

# _pol = Polytope(Fill(HEX_AXIS,3)...)
# _rt_reffe = TensorReferenceFE(_pol,DivConformity(),Float64,0)
# _U = FESpace(model,_rt_reffe) # trial space


u_cf = CellField(u_exact,Ω)
u_int = interpolate(u_exact,U)

n = get_normal_vector(Ω)
n_norm = Operation(n -> n⋅n)(n)
g = evaluate(n/n_norm,get_cell_points(Ω))./1
# u_n = u_cf - (u_cf⋅n)*n/n_norm
u_n = u_cf - (u_cf⋅(n/n_norm))*n/n_norm


# project exact solution
a(u,v) = ∫( u⋅v )dΩ
l(v) = ∫( u_n⋅v  )dΩ
op = AffineFEOperator(a,l,U,V)
u_proj = solve(LUSolver(),op)

u_proj_n = u_proj - (u_proj⋅(n/n_norm))*n/n_norm

writevtk(Ω,joinpath(out_dir,"constants_sphere"),
                cellfields=["u_cf"=>u_cf,"u_int"=>u_int,
                            "u_proj"=>u_proj,
                            "u_n" => u_n,
                            "u_proj_n" => u_proj_n,
                            "err_p" => u_n-u_proj_n,
                            "err_i" => u_int-u_cf],append=false)


l2(u_n-u_proj_n,dΩ)


## spherical velocity
U0 = 1.0
# Initial velocity
function u0(xyz)
  θϕr = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  u = U0#*cos(ϕ)
  spherical_to_cartesian_matrix(θϕr)⋅VectorValue(1.0,2.0,3.0)
end

_u_cf = CellField(u0,Ω)
_u_int = interpolate(u0,U)


# project exact solution
a(u,v) = ∫( u⋅v )dΩ
l(v) = ∫( u0⋅v  )dΩ
op = AffineFEOperator(a,l,U,V)
_u_proj = solve(LUSolver(),op)

writevtk(Ω,joinpath(out_dir,"constants_sphere"),
                cellfields=["u_cf"=>_u_cf,"u_int"=>_u_int,
                            "u_proj"=>_u_proj,
                             ],append=false)



### low level cell field -- https://github.com/gridap/Gridap.jl/blob/3b87738524f6eedd187459b769ce572035009678/src/CellData/CellFields.jl#L110
using Gridap.CellData, Gridap.Fields, Gridap.FESpaces, Gridap.Algebra, Gridap.Geometry
using Gridap.Helpers, Gridap.Arrays, Gridap.Polynomials, Gridap.ReferenceFEs
using FillArrays

trian = Triangulation(model)

s = size(get_cell_map(trian))
_cell_field = Fill(GenericField(u_exact),s)
ucf = GenericCellField(_cell_field,trian,PhysicalDomain())

println( evaluate(ucf,get_cell_points(Ω))./1 )

#### low level H div element
order = 0

et = Float64
p=HEX
D = num_dims(p)
is_n_cube(p)
prebasis = QCurlGradMonomialBasis{D}(et,order)

nf_nodes, nf_moments = Gridap.ReferenceFEs._RT_nodes_and_moments(et,p,order,GenericField(identity))
face_own_dofs = Gridap.ReferenceFEs._face_own_dofs_from_moments(nf_moments)
face_dofs = face_own_dofs

dof_basis = MomentBasedDofBasis(nf_nodes, nf_moments)

ndofs = num_dofs(dof_basis)



#### low level interpolate
object(x) = VectorValue(1.0, 2.0, 3.0 )

free_vals = zero_free_values(U)
s = get_fe_dof_basis(U) ### dof basis from Hdiv
trian = get_triangulation(s)
f = CellField(object,trian,DomainStyle(s)) ## physical vs. ref domain

println( evaluate(f,get_cell_points(Ω))./1  ) ### these are correct

cell_vals = s(f) ### only bit from Hdiv
println(cell_vals./1)

cmaps = get_cell_map(model)

U.fe_basis.cell_basis[1]
evaluate(U.fe_basis.cell_basis[1],VectorValue(0.0,0.0))

s2 = map(cm -> RaviartThomasRefFE(Float64,QUAD,0;phi=cm),cmaps)

cell_vals_2 = map(evaluate,s2,CellData.get_data(f))
println(cell_vals_2./1)


dirichlet_vals = zero_dirichlet_values(U)


cell_dofs = get_cell_dof_ids(U)
cache_vals = array_cache(cell_vals)
cache_dofs = array_cache(cell_dofs)
cells = 1:length(cell_vals)

vals = getindex!(cache_vals,cell_vals,1)
dofs = getindex!(cache_dofs,cell_dofs,1)


  for cell in cells
    vals = getindex!(cache_vals,cell_vals,cell)
    dofs = getindex!(cache_dofs,cell_dofs,cell)
    for (i,dof) in enumerate(dofs)
      val = vals[i]
      if dof > 0
        println("here")
        free_vals[dof] = val
      elseif dof < 0
        dirichlet_vals[-dof] = val
      else
        @unreachable "dof ids either positive or negative, not zero"
      end
    end
  end

  (free_vals,dirichlet_vals)
  free_vals


diri_values = get_dirichlet_dof_values(U)

@check eltype(free_vals) == eltype(diri_values)
cell_dof_ids = get_cell_dof_ids(U)


g = PosNegReindex(free_vals,diri_values)


_cell_vals = lazy_map(Broadcasting(PosNegReindex(free_vals,diri_values)),cell_dof_ids)
cell_field = CellField(U,_cell_vals)

println( evaluate(cell_field,get_cell_points(Ω))./1 )

uh =  SingleFieldFEFunction(cell_field,_cell_vals,free_vals,diri_values,U)

writevtk(Ω,joinpath(out_dir,"constants_sphere"),
                cellfields=["u_cf"=>u_cf,"u_int"=>u_int,
                            "u_proj"=>u_proj,"uh"=>uh],append=false)
