using Gridap
include("../src/initialise.jl")
include("../src/Visualization/VisualizationData.jl")

################################################################################
##### analytical function in ambient space
################################################################################
function tangent_f(vec_func::Function)
  function _tangent_f(X)
    vec = vec_func(X) # VectorValue(-X[2],X[1],0)
    # vec = VectorValue(X[1]^2,X[1]*X[2],X[3])

    normal_vec = 1/sqrt(X[1]*X[1] + X[2]*X[2] + X[3]*X[3])*VectorValue(X[1],X[2],X[3])
    normal_comp = (vec⋅normal_vec)*normal_vec
    tangent_comp = vec - normal_comp

    @assert norm(normal_vec) ≈ 1.0 # check length of normal vector = 1
    @assert dot(normal_comp, tangent_comp) <= 9.9e-15 # check normal and tangent components are perpendicular

    tangent_comp

  end
end

################################################################################
##### Map analytical vector valued function from X,Y,Z -> α,β
##### 1. Map a point α,β -> X,Y,Z
##### 2. Compute u(X,Y,Z)
##### 3. Compute J cooresponding to α,β -> X,Y,Z
##### 4. Evaluate v(α,β) = pinv(J)⋅u(X,Y,Z)
## Note, this is a coordinate transform of an analytic function, so do not need to
## apply the piola transform (I think)
## Note, the inverse map X,Y,Z -> α,β is ill defined (det(J) = 0). This is why
## we need to use the (left) pseudo-inverse of the forward map
################################################################################
function u_ambient_vector(p::Int,uX::Function)
  function _u(αβ)

    cmap = PanelRotationField(r1p_3D[p]) ∘ SigmaField(r) ∘ GnomonicField()
    inv_cmap = InvGnomonicField() ∘ InvSigmaField(r)  ∘ PanelRotationField(rp1_3D[p])

    XYZ = cmap(αβ)

    Jt = ∇(cmap)
    J = Operation(transpose)(Jt)
    Jt_x = Jt(αβ)
    J_x = J(αβ)

    pinvJ = inv(Jt_x⋅J_x)⋅Jt_x

    pinvJ ⋅ uX(XYZ)

  end
end

################################################################################
##### Interoplate cell-wise array of generic fields
################################################################################
function Gridap.FESpaces._cell_vals(fs::SingleFieldFESpace,object::AbstractArray{<:GenericField})
  println("my cell vals")
  s = get_fe_dof_basis(fs)
  trian = get_triangulation(s)
  f = CellData.GenericCellField(object,trian,PhysicalDomain())
  s(f)
end

function Gridap.FESpaces._cell_vals(fs::SingleFieldFESpace,cf_ambient::CellField)
  println("my cell vals")
  s = get_fe_dof_basis(fs)
  s(cf_ambient)
end




using Plots
using LaTeXStrings
# set up plot attributes
markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
_colors = palette(:tab10)
n_ref = [0,1,2,3,4]


### New cubed sphre model: test convergence over series of refined models
manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
models = get_refined_models(manifold_model,4)
vec_func1(X) = VectorValue(-X[2],X[1],0)
vec_func2(X) = VectorValue(-X[1]*X[3],-X[2]*X[3],X[3]*X[3]-r^2)
vec_func3(X) = VectorValue(X[1]*X[1],X[2]*X[2],X[3]*X[3])
vec_func4(X) = VectorValue(X[1]^2,X[1]*X[2],X[3])

vec_funs = [vec_func1, vec_func2, vec_func3, vec_func4]

legendinf = [L"u = (-y,x,0)",L"v = (-xz,-yz,z^2-r^2)",
            L"w = (x^2,y^2,z^2)", L"q = (x^2,xy,z)"]

plot()
for j in 1:length(vec_funs)
  errs = []
  vec = vec_funs[j]
  for i in 1:length(models)
    println(i)
    manifold_model = models[i]

    ## extract domain information
    panel_ids = get_panel_ids(manifold_model)
    ambient_model = get_ambient_model(manifold_model)

    Ω_parametric = Triangulation(manifold_model)
    Ω_ambient = Triangulation(ambient_model)

    ## interpolate analytic function into ambient space
    RT_ambient = FESpace(ambient_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
    analytic_u_ambient = interpolate(tangent_f(vec),RT_ambient)

    ## interpolate mapped analytic function into parametric space
    RT = FESpace(manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
    cell_field = map(p->GenericField(u_ambient_vector(p,tangent_f(vec))),panel_ids)
    uh = interpolate(cell_field,RT)

    ## map parametric FEFunction back to ambient space
    mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids)
    inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

    Jt = lazy_map(Broadcasting(gradient),mapping)
    J = lazy_map(Operation(transpose),Jt)

    _uh = change_domain(uh.cell_field,ReferenceDomain(),PhysicalDomain())

    _cf_mapped = lazy_map(Broadcasting(⋅),J,get_data(_uh))
    cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,inv_mapping)

    cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() ) # ambient cell field

    ### Convert ambient cellfield into a FEFunction by evaluating the field at the dofs
    uh_ambient = interpolate(cf_ambient,RT_ambient)

    writevtk(Ω_ambient,dir*"/ambient_tangent_$i",
          cellfields=["f"=>tangent_f(vec),"u_ambient_cf"=>cf_ambient,"uh_ambient"=>uh_ambient,
          "eh"=>analytic_u_ambient-uh_ambient,
          "u_analytic_ambient"=>analytic_u_ambient],append=false)

    ## compute error
    e = analytic_u_ambient - uh_ambient
    dΩ = Measure(Ω_ambient,2)
    push!(errs, l2(e,dΩ) )
  end
  plot!(n_ref[2:end],errs[2:end],
  lw=4,ms=6,
  c=_colors[j],
  markershape=markers[j],
   label = legendinf[j]
      )

end
plot!(yscale=:log10,framestyle=:box,
  # title = "interpolation error",
  xlabel=L"n",
  ylabel=L"L^2(u_h - u)"
  )
plot!(show=true,legend=:bottomleft)

plot!(xtickfontsize=12,ytickfontsize=12,
legendfontsize=12,guidefontsize=18)

plot!(ylimits=(1e-15,1e-2))
savefig(plotsdir()*"/interpolation_error")
