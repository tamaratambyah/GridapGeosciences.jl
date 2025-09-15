using Gridap.Algebra, Gridap.FESpaces
ls = LUSolver()

### Hdiv
massHdiv(u,v) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ
AHdiv = assemble_matrix(massHdiv,U,V)
nsHdiv = numerical_setup(symbolic_setup(ls,AHdiv),AHdiv)
bHdiv = allocate_in_domain(AHdiv)
assemHdiv = SparseMatrixAssembler(V,V)

### L2
massL2(u,v) = ∫( u*v*meas_cf  )dΩ
AL2 = assemble_matrix(massL2,P,Q)
nsL2 = numerical_setup(symbolic_setup(ls,AL2),AL2)
bL2 = allocate_in_domain(AL2)
assemL2 = SparseMatrixAssembler(Q,Q)

### H1
AH1 = assemble_matrix(massL2,H,R)
bH1 = allocate_in_domain(AH1)
assemAH1 = SparseMatrixAssembler(H,R)
assembH1 = SparseMatrixAssembler(R,R)
nsH1 = numerical_setup(symbolic_setup(ls,AH1),AH1)

### MultiFieldFESpace Hdiv × L2
AX = assemble_matrix(massX,X_prog,Y_prog)
bX = allocate_in_domain(AX)
assemX = SparseMatrixAssembler(Y_prog,Y_prog)
nsX = numerical_setup(symbolic_setup(ls,AX),AX)


### Vorticity functions
perp_matrix_cf = CellField(analytic_perp_matrix,Ω_panel)
biformq((u,p),q,w) = ∫( q*p*w*meas_cf  )dΩ
# liformq((u,p),w) =  ∫( cor_cf*w*meas_cf  )dΩ + ∫( (perp_matrix_cf⋅u)⋅∇(w)  )dΩ
liformq((u,p),w) = ∫( (perp_matrix_cf⋅u)⋅∇(w)  )dΩ
liformq_const(w) =  ∫( cor_cf*w*meas_cf  )dΩ
qconst = assemble_vector(liformq_const,R)

aq(xh) = (q,w) -> biformq(xh,q,w)
bq(xh) = (w) -> liformq(xh,w)

### Mass flux functions
liformF((u,p),v) =  ∫( p*(u⋅(metric_cf⋅v))*meas_cf   )dΩ
bF(xh) = v -> liformF(xh,v)

### Bernoulli functions
# liformΦ((u,p),ψ) =  ∫( gravity*(p+b_cf)*ψ*meas_cf  )dΩ + ∫( 0.5*( u ⋅(metric_cf⋅u) )ψ*meas_cf  )dΩ
liformΦ((u,p),ψ) =  ∫( gravity*p*ψ*meas_cf  )dΩ + ∫( 0.5*( u ⋅(metric_cf⋅u) )ψ*meas_cf  )dΩ
liformΦ_const(ψ) =  ∫( gravity*b_cf*ψ*meas_cf  )dΩ
Φconst =  assemble_vector(liformΦ_const,Q)

bΦ(xh) = ψ -> liformΦ(xh,ψ)


### Prognostic functions
resdepth(F,Φ,q,(v,r)) = ∫( -r*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ
resvelocity((u,p),F,Φ,q,(v,r)) =
                    (  ∫( Φ*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                    - ∫( q*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                    + ∫( τ*(u⋅∇(q))*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
        )

massX((u,p),(v,r)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ +  ∫( (p*r)*meas_cf )dΩ
resX((u,p),F,Φ,q,(v,r)) = resdepth(F,Φ,q,(v,r))  + resvelocity((u,p),F,Φ,q,(v,r))



## work vectors
_F = allocate_in_domain(AHdiv)
q = allocate_in_domain(AH1)
Φ = allocate_in_domain(AL2)
x = allocate_in_domain(AX)

# vorticity solve
function compute_vorticity!(q,bH1,AH1,nsH1,assembH1,assemAH1,xh,H,R,ls)
  # fill!(AH1,0.0)
  # fill!(bH1,0.0)
  fill!(q,0.0)

  # aq(q,w) = biformq(xh,q,w)
  # bq(w) = liformq(xh,w)

  assemble_matrix!(aq(xh),AH1,assemAH1,H,R)
  assemble_vector!(bq(xh),bH1,assembH1,R)
  axpby!(1,qconst,1,bH1)  # add the constant term

  numerical_setup!(nsH1,AH1) # redo numerical set up


  solve!(q,nsH1,bH1)
  return FEFunction(H,q)
end



# mass flux solve
function compute_mass_flux!(_F,bHdiv,AHdiv,nsHdiv,assemHdiv,xh,U,V)
  # fill!(bHdiv,0.0)
  fill!(_F,0.0)
  assemble_vector!(bF(xh),bHdiv,assemHdiv,V)

  solve!(_F,nsHdiv,bHdiv)
  return FEFunction(U,_F)
end




# Bernoulli potential
function compute_bernoulli!(Φ,bL2,AL2,nsL2,assemL2,xh,P,Q)
  # fill!(bL2,0.0)
  fill!(Φ,0.0)
  assemble_vector!(bΦ(xh),bL2,assemL2,Q)
  axpby!(1,Φconst,1,bL2) # add the constant term


  solve!(Φ,nsL2,bL2)
  return FEFunction(P,Φ)
end


function stage1!(x,bX,AX,bXn,nsX,assemX,xn,F,Φ,q,X,Y)
  # fill!(bX,0.0)
  # b1(y) = massX(xn,y) + dt*resX(xn,F,Φ,q,y)
  b1(y) =  dt*resX(xn,F,Φ,q,y)
  assemble_vector!(b1,bX,assemX,Y)

  axpby!(1.0,bXn,1.0,bX)


  fill!(x,0.0)
  solve!(x,nsX,bX)
  return FEFunction(X,x)
end


function stage2!(x,bX,AX,bXn,nsX,assemX,xn,x1,F,Φ,q,X,Y)
  # fill!(bX,0.0)
  # b2(y) = 0.75*massX(xn,y) + 0.25*( massX(x1,y) + dt*resX(x1,F,Φ,q,y) )
  b2(y) =  massX(x1,y) + dt*resX(x1,F,Φ,q,y)
  assemble_vector!(b2,bX,assemX,Y)

  axpby!(0.75,bXn,0.25,bX)

  fill!(x,0.0)
  solve!(x,nsX,bX)
  return FEFunction(X,x)
end


function stage3!(x,bX,AX,bXn,nsX,assemX,xn,x2,F,Φ,q,X,Y)
  # fill!(bX,0.0)
  # b3(y) = (1/3)*massX(xn,y) + (2/3)*( massX(x2,y) + dt*resX(x2,F,Φ,q,y) )
  b3(y) =  massX(x2,y) + dt*resX(x2,F,Φ,q,y)
  assemble_vector!(b3,bX,assemX,Y)

  axpby!((1/3),bXn,(2/3),bX)

  fill!(x,0.0)
  solve!(x,nsX,bX)
  return FEFunction(X,x)
end

dir = datadir("Transient_test")
!isdir(dir) && mkdir(dir)

xh0 = interpolate_everywhere([u_contra_h,h_h],X_prog(0.0))

qh = compute_vorticity!(q,bH1,AH1,nsH1,assembH1,assemAH1,xh0,H,R,ls)
Fh = compute_mass_flux!(_F,bHdiv,AHdiv,nsHdiv,assemHdiv,xh0,U,V)
Φh = compute_bernoulli!(Φ,bL2,AL2,nsL2,assemL2,xh0,P,Q)

Enstropys = Float64[]
Energys = Float64[]
Masss = Float64[]

ens0 = sum(∫( (qh*qh*xh0[2])*meas_cf  )dΩ)
energy0 = sum(∫( (0.5*xh0[2]*( xh0[1] ⋅(metric_cf⋅xh0[1])) + 0.5*gravity*xh0[2]*xh0[2] )*meas_cf )dΩ)
mass0 = sum( ∫( xh0[2]*meas_cf )dΩ  )
push!(Enstropys,ens0)
push!(Energys,energy0)
push!(Masss,mass0)

# bn(y) = massX(xh0,y)
bn(x) = y -> massX(x,y)
bXn = assemble_vector(bn(xh0),assemX,Y_prog)


cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

vort = qh*xh0[2] - cor_cf
labels = ["q","F","Phi","u","p","vort"]
panel_cfs = [qh,Fh,Φh, covarient_basis_cf⋅xh0[1], xh0[2],vort]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)


for nsteps in collect(1:1000)
  t = nsteps*dt
  println(t)

  # stage 1
  # qh = compute_vorticity!(q,bH1,AH1,nsH1,assembH1,assemAH1,xh0,H,R,ls)
  # Fh = compute_mass_flux!(_F,bHdiv,AHdiv,nsHdiv,assemHdiv,xh0,U,V)
  # Φh = compute_bernoulli!(Φ,bL2,AL2,nsL2,assemL2,xh0,P,Q)
  xh1 = stage1!(x,bX,AX,bXn,nsX,assemX,xh0,Fh,Φh,qh,X_prog,Y_prog)

  # stage 2
  qh = compute_vorticity!(q,bH1,AH1,nsH1,assembH1,assemAH1,xh1,H,R,ls)
  Fh = compute_mass_flux!(_F,bHdiv,AHdiv,nsHdiv,assemHdiv,xh1,U,V)
  Φh = compute_bernoulli!(Φ,bL2,AL2,nsL2,assemL2,xh1,P,Q)
  xh2 = stage2!(x,bX,AX,bXn,nsX,assemX,xh0,xh1,Fh,Φh,qh,X_prog,Y_prog)

  # stage 3
  qh = compute_vorticity!(q,bH1,AH1,nsH1,assembH1,assemAH1,xh2,H,R,ls)
  Fh = compute_mass_flux!(_F,bHdiv,AHdiv,nsHdiv,assemHdiv,xh2,U,V)
  Φh = compute_bernoulli!(Φ,bL2,AL2,nsL2,assemL2,xh2,P,Q)
  xh3 = stage2!(x,bX,AX,bXn,nsX,assemX,xh0,xh2,Fh,Φh,qh,X_prog,Y_prog)

  # final diagnostics
  qh = compute_vorticity!(q,bH1,AH1,nsH1,assembH1,assemAH1,xh3,H,R,ls)
  Fh = compute_mass_flux!(_F,bHdiv,AHdiv,nsHdiv,assemHdiv,xh3,U,V)
  Φh = compute_bernoulli!(Φ,bL2,AL2,nsL2,assemL2,xh3,P,Q)

  # println("casimirs")

  ens = sum(∫( (qh*qh*xh3[2])*meas_cf  )dΩ)
  energy = sum(∫( (0.5*xh3[2]*( xh3[1] ⋅(metric_cf⋅xh3[1])) + 0.5*gravity*xh3[2]*xh3[2] )*meas_cf )dΩ)
  mass = sum( ∫( xh3[2]*meas_cf )dΩ  )

  push!(Enstropys,ens)
  push!(Energys,energy)
  push!(Masss,mass)

  println("enstropy:  ", abs(ens-ens0)/ens0)
  println("energy:  ", abs(energy-energy0)/energy0)
  println("mass:  ", abs(mass-mass0)/mass0)

  vort = qh*xh3[2] - cor_cf

  labels = ["q","F","Phi","u","p","vort"]
  panel_cfs = [qh,Fh,Φh, covarient_basis_cf⋅xh3[1], xh3[2],vort]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)
  writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)


  xh0 =  FEFunction(X_prog,x)
  assemble_vector!(bn(xh0),bXn,assemX,Y_prog)

end

make_pvd(dir,"solT",1)


ts = dt*collect(0:length(Masss)-1)

ms_rel = abs.(Masss.-Masss[1])./Masss[1]
Es_rel = abs.(Energys.-Energys[1])./Energys[1]
Enst_rel = abs.(Enstropys.-Enstropys[1])./Enstropys[1]

plot()
plot!(ts[2:end],ms_rel[2:end],lw=3,label="mass")
plot!(yaxis=:log,xlabel="t",ylabel=L"|x_t-x_0|/x_0")

plot()
plot!(ts[2:end],Es_rel[2:end],lw=3,label="energy")
plot!(yaxis=:log,xlabel="t",ylabel=L"|x_t-x_0|/x_0")
# savefig(plotsdir()*"/wave_transient_conservation")

plot()
plot!(ts[2:end],Enst_rel[2:end],lw=3,label="enstropy")
plot!(yaxis=:log,xlabel="t",ylabel=L"|x_t-x_0|/x_0")
