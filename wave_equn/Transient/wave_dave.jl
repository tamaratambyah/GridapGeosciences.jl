
using Gridap
using GridapSolvers
using GridapGeosciences
using LinearAlgebra
using DrWatson

include("createpvd.jl")

# Initial depth
function h₀(xyz)
  x,y,z = xyz
  xᵢ,yᵢ,zᵢ = (1,0,0)
  1.0 + 0.0001*exp(-5*((xᵢ-x)^2+(yᵢ-y)^2+(zᵢ-z)^2))
end

# Initial velocity
u₀(xyz) = zero(xyz)

function new_vtk_step(Ω,file,hn,un)
  createvtk(Ω,
            file,
            cellfields=["hn"=>hn, "un"=>un],
            append=false)
end

function _ssrk2_update!(res,MM,MMchol,DIV,α,a,b)
  mul!(res, MM, a)            # res <- MM*a
  mul!(res, DIV, b, α, 1.0)   # res <- 1.0*res + α*DIV*b
  solve!(res,MMchol,res)           # res <- inv(MM)*res
end


"""
  Solves the wave equation using a 2nd order
  Strong-Stability-Preserving Runge-Kutta
  explicit time integration scheme (SSRK2)

  g : acceleration due to gravity
  H : (constant) reference layer depth
  T : [0,T] simulation interval
  N : number of time subintervals
"""
function solve_wave_equation_ssrk2(
      model,order,degree,g,H,T,N;
      mass_matrix_solver::Gridap.Algebra.LinearSolver=Gridap.Algebra.BackslashSolver(),
      write_results=false,
      out_dir=datadir("wave_eq_ncells_$(num_cells(model))_order_$(order)_ssrk2"),
      out_period=N/10)

  !isdir(out_dir) && mkdir(out_dir)
  pvd = createpvd(out_dir)

  RT=ReferenceFE(raviart_thomas,Float64,order)
  DG=ReferenceFE(lagrangian,Float64,order)
  V = FESpace(model,RT; conformity=:Hdiv)
  Q = FESpace(model,DG; conformity=:L2)

  U = TrialFESpace(V)
  P = TrialFESpace(Q)

  Ω  = Triangulation(model)
  dΩ = Measure(Ω,degree)
  dω = Measure(Ω,degree,ReferenceDomain())

  # Build mass matrix in RT and DG spaces
  # and their sparse Cholesky factors
  amm(u,v) = ∫(v⋅u)dΩ
  RTMM=assemble_matrix(amm,U,V)
  L2MM=assemble_matrix(amm,P,Q)
  RTMMchol=numerical_setup(symbolic_setup(mass_matrix_solver,RTMM),RTMM)
  L2MMchol=numerical_setup(symbolic_setup(mass_matrix_solver,L2MM),L2MM)

  # Build g*DIV(v)*h and H*q*DIV(u)
  ad(h,v) = ∫(DIV(v)*h)dω
  divvh=assemble_matrix(ad,P,V)
  adt(u,q) = ∫(q*DIV(u))dω
  qdivu=assemble_matrix(adt,U,Q)

  # Interpolate initial condition into FE spaces
  hn=interpolate_everywhere(h₀,P); hnv=get_free_dof_values(hn)
  un=interpolate_everywhere(u₀,U); unv=get_free_dof_values(un)

  pvd[0] = new_vtk_step(Ω,joinpath(out_dir,"n_t0"),hn,un)

  # Allocate work space vectors
  h1v = similar(get_free_dof_values(hn))
  h2v = similar(h1v)
  u1v = similar(get_free_dof_values(un))
  u2v = similar(u1v)

  dt  = T/N
  dtg = dt*g
  dtH = dt*H



  @time for step=1:N
      # 1st step
      # inv(L2MM)*(L2MM*hnv - dt*Hqdivu*unv)
      # inv(RTMM)*(RTMM*unv + dt*gdivvh*hnv)
      _ssrk2_update!(h1v, L2MM, L2MMchol, qdivu,-dtH, hnv, unv)
      _ssrk2_update!(u1v, RTMM, RTMMchol, divvh, dtg, unv, hnv)

      # 2nd step
      # inv(L2MM)*(L2MM*h1v - dt*Hqdivu*u1v)
      # inv(RTMM)*(RTMM*u1v + dt*gdivvh*h1v)
      _ssrk2_update!(h2v, L2MM, L2MMchol, qdivu, -dtH, h1v, u1v)
      _ssrk2_update!(u2v, RTMM, RTMMchol, divvh,  dtg, u1v, h1v)

      # Averaging steps
      hnv .= 0.5 .* ( hnv .+ h2v )
      unv .= 0.5 .* ( unv .+ u2v )



    println(step)
    pvd[step] = new_vtk_step(Ω,joinpath(out_dir,"n_t$(step)"),hn,un)


  end
  # end

  pvd[N] = new_vtk_step(Ω,joinpath(out_dir,"n_t$(N)"),hn,un)

  un,hn,out_dir


end

model=CubedSphereDiscreteModel(4)
g=1.0
H=1.0
T=π
N=2000
order=0
degree=6
un,hn,out_dir =
  solve_wave_equation_ssrk2(model,order,degree,g,H,T,N;write_results=true,out_period=10)

make_pvd(out_dir,"n__","dave_sol",1)
