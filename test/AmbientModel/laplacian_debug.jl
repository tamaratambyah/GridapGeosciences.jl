
using Gridap
using GridapGeosciences
using Test

dir = @__DIR__


radius = 1
n_ref_lvls = 1
ambient_model = CubedSphereAmbientDiscreteModel(
  radius; num_initial_uniform_refinements=n_ref_lvls)
panel_model = ambient_model.panel_model


################################################################################
#### f = -sin(ϕ)
#### This is Williamson2, equation (92), with u₀=0, α = 0 (not panel coord)
#### The analytic laplacian is Williamson2, equation (94), with u₀=0, Ω=0, α = 0 (not panel coord)
################################################################################
function analytic_f_θϕr(θϕr)
  θ,ϕ,r = θϕr
  -sin(ϕ)
end

function analytic_lap_θϕr(θϕr)
  θ,ϕ,r = θϕr
  2*sin(ϕ)
end

### Ambient model
function ambient_fX(xyz)
  θϕr   = xyz2θϕr(xyz)
  analytic_f_θϕr(θϕr)
end
ambient_lapX(x) = (∇⋅ ∇(ambient_fX))(x)
Ω_ambient = Triangulation(ambient_model)
ambient_fX_cf = CellField(ambient_fX,Ω_ambient)
ambient_lapX_cf = CellField(ambient_lapX,Ω_ambient)

writevtk(Ω_ambient,dir*"/ambient_trian",
cellfields=["f"=>ambient_fX_cf,"slap"=>ambient_lapX_cf,"slap_analytic"=>ambient_lapX,
            "e_slap"=>ambient_lapX_cf-ambient_lapX ],
append=false)


### Parametric model
function panel_fX(forward_map)
  function _f(α)
    xyz = forward_map(α)
    θϕr   = xyz2θϕr(xyz)
    analytic_f_θϕr(θϕr)
  end
end

function panel_surflap_analytic(forward_map)
  function _f(α)
    xyz = forward_map(α)
    θϕr   = xyz2θϕr(xyz)
    analytic_lap_θϕr(θϕr)
  end
end


Ω_panel = Triangulation(panel_model)
panel_f_cf = ParametricCellField(panel_fX,Ω_panel)
panel_slap_cf =  ParametricCellField(surflap(panel_fX),Ω_panel)
analytic_slap_cf = ParametricCellField(panel_surflap_analytic,Ω_panel)

writevtk_with_cell_geomap(geo_map_func(Ω_panel),Ω_panel,dir*"/panel_trian",
cellfields=["f"=>panel_f_cf,"slap"=>panel_slap_cf,"slap_analytic"=>analytic_slap_cf,
            "e_slap"=>panel_slap_cf-analytic_slap_cf ],
append=false)


################################################################################
#### f = xyz
#### No ananlytic solution for laplacian
################################################################################

## Ambient model
function ambient_fX(xyz)
  x,y,z = xyz
  x*y*z
end
ambient_lapX(x) = (∇⋅ ∇(ambient_fX))(x)
Ω_ambient = Triangulation(ambient_model)
ambient_fX_cf = CellField(ambient_fX,Ω_ambient)
ambient_lapX_cf = CellField(ambient_lapX,Ω_ambient)

writevtk(Ω_ambient,dir*"/ambient_trian",
cellfields=["f"=>ambient_fX_cf,"slap"=>ambient_lapX_cf ],
append=false)


## Panel model
function panel_fX(forward_map)
  function _f(αβ)
    x = forward_map(αβ)
    x[1]*x[2]*x[3]
  end
end

Ω_panel = Triangulation(panel_model)
panel_f_cf = ParametricCellField(panel_fX,Ω_panel)
panel_slap_cf =  ParametricCellField(surflap(panel_fX),Ω_panel)

writevtk_with_cell_geomap(geo_map_func(Ω_panel),Ω_panel,dir*"/panel_trian",
cellfields=["f"=>panel_f_cf,"slap"=>panel_slap_cf ],
append=false)
