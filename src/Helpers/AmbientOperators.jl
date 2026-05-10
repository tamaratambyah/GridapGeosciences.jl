function panel_f(ambient_f::Function,forward_map::Field)
  function _f(αβ)
    x = forward_map(αβ)
    ambient_f(x)
  end
end
panel_f(f::Function) = m -> panel_f(f,m)

### How to ensure the function is in the tangent space?
### We cannot remove the normal component, since this is not general to 3D
### Rely on an informed using the ambient_vec ∈ TₚS
function panel_vec(ambient_vec::Function,forward_map::Field)
  function _v(αβ)
    x = forward_map(αβ)
    forward_pinv_jacobian(forward_map,αβ) ⋅ ambient_vec(x)
  end
end
panel_vec(f::Function) = m -> panel_vec(f,m)

ambient_sgrad(f::Function) = m -> ambient_sgrad(f,m)
function ambient_sgrad(f::Function,inverse_map::Field)
  function _gradf(x)
    αβ = inverse_map(x)
    m = ForwardMap(inverse_map.panel,inverse_map.radius)
    sgrad(panel_f(f),m)(αβ)
  end
end

ambient_surflap(f::Function) = m -> ambient_surflap(f,m)
function ambient_surflap(f::Function,inverse_map::Field)
  function _lapf(x)
    αβ = inverse_map(x)
    m = ForwardMap(inverse_map.panel,inverse_map.radius)
    surflap(panel_f(f),m)(αβ)
  end
end

ambient_surfdiv(f::Function) = m -> ambient_surfdiv(f,m)
function ambient_surfdiv(f::Function,inverse_map::Field)
  function _sdivf(x)
    αβ = inverse_map(x)
    m = ForwardMap(inverse_map.panel,inverse_map.radius)
    surfdiv(panel_vec(f),m)(αβ)
  end
end
