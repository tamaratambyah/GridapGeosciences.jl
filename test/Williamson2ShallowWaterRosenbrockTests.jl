module Williamson2ShallowWaterRosenbrockTests

using FillArrays
using Test
using WriteVTK
using Gridap
using GridapGeosciences

# Solves the steady state Williamson2 test case for the shallow water equations on a sphere
# of physical radius 6371220m. Involves a modified coriolis term that exactly balances
# the potential gradient term to achieve a steady state
# reference:
# D. L. Williamson, J. B. Drake, J. J.HackRüdiger Jakob, P. N.Swarztrauber, (1992)
# J Comp. Phys. 102 211-224

# Constants of the Williamson2 test case
const α  = π/4.0              # deviation of the coriolis term from zonal forcing
const U₀ = 38.61068276698372  # velocity scale
const H₀ = 2998.1154702758267 # mean fluid depth

# Modified coriolis term
function f₀(xyz)
   θϕr   = xyz2θϕr(xyz)
   θ,ϕ,r = θϕr
   2.0*Ωₑ*( -cos(θ)*cos(ϕ)*sin(α) + sin(ϕ)*cos(α) )
end

# Initial velocity
function u₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  u     = U₀*(cos(ϕ)*cos(α) + cos(θ)*sin(ϕ)*sin(α))
  v     = -U₀*sin(θ)*sin(α)
  spherical_to_cartesian_matrix(θϕr)⋅VectorValue(u,v,0)
end

# Initial fluid depth
function h₀(xyz)
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  h  = -cos(θ)*cos(ϕ)*sin(α) + sin(ϕ)*cos(α)
  H₀ - (rₑ*Ωₑ*U₀ + 0.5*U₀*U₀)*h*h/g
end

l2_err_u = [0.012427571457879082, 0.0029686291912212367, 0.0007555653294017498]
l2_err_h = [0.005643644188108772, 0.0014535790630939244, 0.0003672222114589216]

order  = 1
degree = 4

λ = 0.5 # magnitude of the descent direction of the implicit solve (neutrally stable for 0.5)

for i in 1:3
  n      = 2*2^i
  nstep  = 5*n
  Uc     = sqrt(g*H₀)
  dx     = 2.0*π*rₑ/(4*n)
  dt     = 0.25*dx/Uc
  println("timestep: ", dt)   # gravity wave time step

  model = CubedSphereDiscreteModel(n; radius=rₑ)
  hf, uf = shallow_water_rosenbrock_time_stepper(model, order, degree,
                                                h₀, u₀, f₀, g, H₀,
                                                λ, dt, 0.0, nstep;
                                                leap_frog=true,
                                                write_solution=false,
                                                write_solution_freq=5,
                                                write_diagnostics=true,
                                                write_diagnostics_freq=1,
                                                dump_diagnostics_on_screen=true)

  Ω     = Triangulation(model)
  dΩ    = Measure(Ω, degree)
  hc    = CellField(h₀, Ω)
  e     = h₀-hf
  err_h = sqrt(sum(∫(e⋅e)*dΩ))/sqrt(sum(∫(hc⋅hc)*dΩ))
  uc    = CellField(u₀, Ω)
  e     = u₀-uf
  err_u = sqrt(sum(∫(e⋅e)*dΩ))/sqrt(sum(∫(uc⋅uc)*dΩ))
  println("n=", n, ",\terr_u: ", err_u, ",\terr_h: ", err_h)

  @test abs(err_u - l2_err_u[i]) < 10.0^-12
  @test abs(err_h - l2_err_h[i]) < 10.0^-12
end

end
