"""
analytic functions defined in terms of the panel coordinates
"""

### f = sin(ϕ) = Z
function f_sin(p)
  function _f(αβ)
    α,β = αβ
    if p == 1 || p == 3
      return RADIUS/rho(αβ)*tan(β)
    elseif p == 2
      return RADIUS/rho(αβ)
    elseif p == 4 || p == 6
      return RADIUS/rho(αβ)*tan(α)
    elseif p == 5
      return -RADIUS/rho(αβ)
    end
  end
end

### f = XYZ
function f_XYZ(p)
  function _f(αβ)
    α,β = αβ
    if p == 1 || p == 5 || p == 6
      return RADIUS^3/rho3(αβ)*tan(α)*tan(β)
    else
      return -RADIUS^3/rho3(αβ)*tan(α)*tan(β)
    end
  end
end

# analytic functions defined on the chart
analytic_funcs = Dict{Symbol,Any}()
analytic_funcs[:sin] = f_sin
analytic_funcs[:XYZ] = f_XYZ

"""
mapped functions defined in terms cartesian coordinates and latlon coordinates
"""

### f = sin(ϕ)
fθϕ(θϕ) = sin(θϕ[2])
fX(XYZ::VectorValue{3}) = XYZ[1]*XYZ[2]*XYZ[3]

# mapped functions from cartesian or latlon
mapped_funcs = Dict{Symbol,Any}()
mapped_funcs[:sin] = panel_to_latlon(fθϕ)
mapped_funcs[:XYZ] = panel_to_cartesian(fX)

"""
Williamson2 stream function, defined using latlon coordinates
"""

fWilliamson(ζ) = θϕ -> - sin(θϕ[2])*cos(ζ) + cos(θϕ[1])*cos(θϕ[2])*sin(ζ)

williamson_funcs = Dict{Symbol,Any}()
williamson_funcs[:z1] = panel_to_latlon(fWilliamson(0))
williamson_funcs[:z2] = panel_to_latlon(fWilliamson(0.05))
williamson_funcs[:z3] = panel_to_latlon(fWilliamson(π/2-0.05))
williamson_funcs[:z4] = panel_to_latlon(fWilliamson(π/2))
