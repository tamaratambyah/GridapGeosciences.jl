include("../AmbientLaplaceBeltrami.jl")

## Serial model: 2D
n_ref_lvls = 4
radius = 1.0
models = get_ambient_refined_models(n_ref_lvls,radius)
AmbientLaplaceBeltrami.main(models)



# ### I do not like having this here, but need to think of a better way
# to compare one error result to the intrinsic approach
ambient_model = models[1]
dir = @__DIR__
p_fe = 2
e_ambient, = AmbientLaplaceBeltrami.laplace_beltrami_solver(
              ambient_model,p_fe,dir,
              AmbientLaplaceBeltrami.fX)


include("../../Laplacian/LaplaceBeltrami.jl")
function fX_ambient(forward_map)
  function _f(αβ)
    x = forward_map(αβ)
    AmbientLaplaceBeltrami.fX(x)
  end
end

panel_model = get_parametric_model(ambient_model)
e_panel, = LaplaceBeltramiTests.laplace_beltrami_solver(
              panel_model,p_fe,dir,
              fX_ambient)

e_comparison = e_ambient - e_panel
println(e_comparison)
@test e_comparison < 1e-12
