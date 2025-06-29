using DrWatson
using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots
include("../../../src/initialise.jl")

function pθϕ_scalar(θϕ)
  θ,ϕ = θϕ
  1 - (sin(ϕ))^2
end

function uθϕ_vector(θϕ)
  θ,ϕ = θϕ
  VectorValue(cos(ϕ),0.0)
end

function uX_vector(X)
  VectorValue(X[1],0.0,0.0)
end

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

p_fe = 1
degree = 2*(p_fe+1)

################################################################################
##### ambient space
################################################################################
ambient_model = get_ambient_model(manifold_model)
Ω_amb = Triangulation(ambient_model)
dΩ_amb = Measure(Ω_amb,degree)

cell_field = map(p->GenericField(u_scalar_latlon2ambient(p,pθϕ_scalar)),panel_ids)
pcf_ambient = CellData.GenericCellField(cell_field,Ω_amb,PhysicalDomain())

# cell_field = map(p->GenericField(u_vector_latlon2ambient(uθϕ_vector)),panel_ids)
# ucf_ambient = CellData.GenericCellField(cell_field,Ω_amb,PhysicalDomain())
ucf_ambient = CellField(uX_vector,Ω_amb)

writevtk(Ω_amb,dir*"/wave_ambient",cellfields=["u"=>ucf_ambient,"p"=>pcf_ambient],append=false)


function u_vector_latlon2parametric(p::Int,uθϕ::Function)
  function _u(αβ)
    function uX(XYZ)
      θϕ = InvSigmaField(RADIUS)(XYZ)

      Jt = gradient(SigmaField(RADIUS)) # 2 x 3
      J = Operation(transpose)(Jt) # 3 x 2

      J(θϕ) ⋅ uθϕ(θϕ)
    end

    cmap = PanelRotationField(r1p_3D[p]) ∘ SigmaField(RADIUS) ∘ GnomonicField()

    XYZ = cmap(αβ)
    θϕ = InvSigmaField(RADIUS)(XYZ)

    Jt = ∇(cmap)
    J = Operation(transpose)(Jt) # 3 x 2
    pinvJ = Operation(pinv)(J) # 2 x 3

    Jt_x = Jt(αβ)
    J_x = J(αβ)
    pinvJ_x = pinvJ(αβ)

    # (2x3) (3x1) = (2x1) ∈ α,β
    return pinvJ_x ⋅ uX(XYZ)
  end
end



################################################################################
##### parametric space
################################################################################
####### check compatibility in parametric space
Ω = Triangulation(manifold_model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω,degree)
dΩg =  Measure(m,Ω,degree)

cell_field = map(p->GenericField(u_scalar_latlon2parametric(p,pθϕ_scalar)),panel_ids)
pcf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

# cell_field = map(p->GenericField(u_vector_latlon2parametric(p,uθϕ_vector)),panel_ids)
cell_field = map(p->GenericField(u_vector_ambient2parametric(p,uX_vector)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

for pid in collect(1:6)
  mask = (panel_ids.== pid)
  Ωp = Triangulation(manifold_model,mask)

  up = change_domain(ucf,Ωp,DomainStyle(ucf))
  pp = change_domain(pcf,Ωp,DomainStyle(pcf))

  writevtk(Ωp,dir*"/wave_p$(pid)",cellfields=["u"=>up,"p"=>pp],append=false)
end


u0 = ucf + surface_gradient(pcf,m)
p0 = pcf + surface_divergence(ucf,m)

### FE problem
V = TestFESpace(Ω,ReferenceFE(raviart_thomas,Float64,p_fe),conformity=:HDiv)
U = TrialFESpace(V)

Q = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,p_fe),conformity=:L2)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V,Q])

biform1(u,v) = ( ∫( u⋅v )dΩg )
biform2(p,q) = ∫( p*q )dΩg
biform3(p,v) =  ∫( -1.0*p*(wave_divergence(v,m))   )dΩ
biform4(u,q) =∫( (surface_divergence(u,m))*q )dΩg

A1 = assemble_matrix(biform1,U,V)
A2 = assemble_matrix(biform2,P,Q)
A3 = assemble_matrix(biform3,P,V)
A4 = assemble_matrix(biform4,U,Q)

liform1(v) = ∫( u0⋅v )dΩg
liform2(q) = ∫(  p0*q  )dΩg
b1 = assemble_vector(liform1,V)
b2 = assemble_vector(liform2,Q)
b = assemble_vector(wave_liform,Y)


wave_biform1((u,p),(v,q)) = ( ∫( u⋅v )dΩg   + ∫( -1.0*p*(wave_divergence(v,m))   )dΩ
                    + ∫( p*q )dΩg + ∫( (surface_divergence(u,m))*q )dΩg
                    )
A = assemble_matrix(wave_biform1,X,Y)
x = A\b
xh = FEFunction(X,x)
uh,ph = xh

wave_biform((u,p),(v,q)) = ( ∫( u⋅v )dΩg
                    + ∫( -1.0*p*(wave_divergence(v,m))   )dΩ
                    + ∫( p*q )dΩg
                    + ∫( (surface_divergence(u,m))*q )dΩg
                    )
wave_liform((v,q)) = ∫( u0⋅v + p0*q  )dΩg

op = AffineFEOperator(wave_biform,wave_liform,X,Y)
uh, ph = solve(LUSolver(),op)


for pid in collect(1:6)
  mask = (panel_ids.== pid)
  Ωp = Triangulation(manifold_model,mask)

  up = change_domain(ucf,Ωp,DomainStyle(ucf))
  pp = change_domain(pcf,Ωp,DomainStyle(pcf))

  uhp = change_domain(uh,Ωp,DomainStyle(uh))
  php = change_domain(ph,Ωp,DomainStyle(ph))
  writevtk(Ωp,dir*"/wave_uhp$(pid)",cellfields=["u"=>up,"uh"=>uhp,"p"=>pp,"ph"=>php],append=false)
end

## map parametric FEFunction back to ambient space
uh_ambient = parametric_cf_2_ambient_vector(manifold_model,uh.cell_field)
ph_ambient = parametric_cf_2_ambient(manifold_model,ph.cell_field)

#### Compute errors
eu =  sum(∫((uh-ucf)⊙(uh-ucf))dΩ)
ep = sum(∫((ph-pcf)⊙(ph-pcf))dΩ)

eu_g =  sum(∫((uh-ucf)⊙(uh-ucf))dΩg)
ep_g = sum(∫((ph-pcf)⊙(ph-pcf))dΩg)

eu_amb = l2(uh_ambient-ucf_ambient,dΩ_amb)
ep_amb = l2(ph_ambient-ucf_ambient,dΩ_amb)
writevtk(Ω_amb,dir*"/wave_ambient",cellfields=["u"=>ucf_ambient,"p"=>pcf_ambient,
"uh"=>uh_ambient,"ph"=>ph_ambient,
"eu"=>ucf_ambient-uh_ambient,"ep"=>pcf_ambient-ph_ambient],append=false)
