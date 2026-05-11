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
function ambient_sgrad(f::Function,forward_map::Field)
  function _gradf(x)
    inverse_map  = InverseMap(forward_map)
    αβ = inverse_map(x)
    sgrad(panel_f(f),forward_map)(αβ)
  end
end

ambient_surflap(f::Function) = m -> ambient_surflap(f,m)
function ambient_surflap(f::Function,forward_map::Field)
  function _lapf(x)
    inverse_map  = InverseMap(forward_map)
    αβ = inverse_map(x)
    surflap(panel_f(f),forward_map)(αβ)
  end
end

ambient_surfdiv(f::Function) = m -> ambient_surfdiv(f,m)
function ambient_surfdiv(f::Function,forward_map::Field)
  function _sdivf(x)
    inverse_map  = InverseMap(forward_map)
    αβ = inverse_map(x)
    surfdiv(panel_vec(f),forward_map)(αβ)
  end
end
