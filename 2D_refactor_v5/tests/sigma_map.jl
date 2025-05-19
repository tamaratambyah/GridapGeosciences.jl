using Gridap
include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)

ambient_model = get_ambient_model(manifold_model)
panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)

quad_parametric = CellQuadrature(Ω_parametric,2)
quad_ambient = CellQuadrature(Ω_ambient,2)

####
# function Fields.pinvJt(Jt::MultiValue{Tuple{D1,D2}}) where {D1,D2}
#   # @check D1 < D2
#   if D1 < D2
#     J = transpose(Jt)
#     return transpose(inv(Jt⋅J)⋅Jt)
#   else
#     println("my inverse")
#     J = transpose(Jt)
#     return transpose(Jt⋅inv(J⋅Jt))
#   end
# end

################################################################################
#### consider sigma map -> only defined on quad points
################################################################################
quad_pt = quad_parametric.cell_point[1]
Jt = ∇(SigmaField(r))
Jt(quad_pt)


#### consider a FEfunction defined in parametric space:
alpha(x) = x[1]
H1 = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,1),conformity=:H1)
uh = interpolate_everywhere(alpha,H1) # 2D ref -> 2D alpha,beta
uh_cf = get_data(uh)
cvals_parametric = lazy_map(evaluate,uh_cf,quad_parametric.cell_point)


m = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids) # 2D alpha,beta -> 3D X,Y,Z sphere

uh_cf_mapped = lazy_map(Broadcasting(∘),uh_cf,m)
# cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,DomainStyle(uh) )
cvals_ambient =  lazy_map(evaluate,uh_cf_mapped,quad_ambient.cell_point)

Jt = lazy_map(Broadcasting(∇),m)
det_Jt = lazy_map(Operation(meas),(Jt))
det_Jtx = lazy_map(evaluate,det_Jt,quad_parametric.cell_point)



























################################################################################
#### debugging
# _pinvJt = lazy_map(Broadcasting(pinvJt),Jt)
pinvJtx_p2 = lazy_map(evaluate,_pinvJt,quad_parametric.cell_point)


#### consider sigma inverse map -> only defined on quad points
c = get_cell_map(ambient_model) # 2D ref -> 3D X,Y,Z
invSigma  = map(x-> InvSigmaField(r), panel_ids) # 3D X,Y,Z -> 2D latlon
C = lazy_map(∘,invSigma,c) # 2D ref -> 2D latlon

invJt = lazy_map(Broadcasting(∇),C)
inv_pinvJt = lazy_map(Broadcasting(pinvJt),invJt)
inv_pinvJtx_p2 =  lazy_map(evaluate,inv_pinvJt,quad_ambient.cell_point)[panel_ids.==2]





m = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(r) ∘ GnomonicField() , panel_ids) # 2D alpha,beta -> 3D X,Y,Z sphere
_uh_mapped = lazy_map(Broadcasting(∘),uh_cf,m)
_uh_cf = CellData.GenericCellField(_uh_mapped,Ω_ambient,DomainStyle(uh) )

cvals_ambient = lazy_map(evaluate,_uh_cf,quad_parametric.cell_point)

cmap = get_cell_map(manifold_model) # 2D ref -> 2D alpha,beta

M = lazy_map(∘,m,cmap) # 2D ref -> 3D X,Y,Z
