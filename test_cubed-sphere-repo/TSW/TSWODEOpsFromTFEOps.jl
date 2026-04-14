

#######################
# PoissonODEOpFromTFEOpCache #
#######################
"""
    struct PoissonODEOpFromTFEOpCache <: GridapType

Structure that stores the `TransientFESpace` and cache of a
`TransientFEOperator`, as well as the jacobian matrices and residual if they
are constant.
"""
mutable struct TSWODEOpFromTFEOpCache <: GridapType
  Us # u, h, B
  Uts
  tfeopcache
  Ys # F, Φ, q
  Yts
  diagnostics
  diagslvrcache
  Bs # b
  Bts
  diagnosticsb
  diagslvrcacheb
  Ws # bhat
  Wts
  diagnosticsw
  diagslvrcachew
  Zs # btilde
  Zts
  diagnosticsz
  diagslvrcachez
  Ts
  Tts
  Bbars
  Bbarts
end

##################
# ODEOpFromTFEOp #
##################
"""
    struct ODEOpFromTFEOp <: ODEOperator end

Wrapper that transforms a `TransientFEOperator` into an `ODEOperator`, i.e.
takes `residual(t, uh, ∂t[uh], ..., ∂t^N[uh], vh)` and returns
`residual(t, us)`, where `us[k] = ∂t^k[us]` and `uf` represents the free values
of `uh`.
"""
struct TSWODEOpFromTFEOp{T} <: TSWODEOperator{T}
  tfeop::TransientTSWOperator
  function TSWODEOpFromTFEOp(tfeop::TransientTSWOperator)
    new{NonlinearODE}(tfeop)
  end
end

# ODEOperator interface
function Gridap.Polynomials.get_order(odeop::TSWODEOpFromTFEOp)
  Gridap.ODEs.get_order(odeop.tfeop)
end

function get_progop(odeop::TSWODEOpFromTFEOp)
  get_progop(odeop.tfeop)
end

function get_diagop(odeop::TSWODEOpFromTFEOp)
  get_diagop(odeop.tfeop)
end

function get_diagopb(odeop::TSWODEOpFromTFEOp)
  get_diagopb(odeop.tfeop)
end

function get_diagopw(odeop::TSWODEOpFromTFEOp)
  get_diagopw(odeop.tfeop)
end

function get_diagopz(odeop::TSWODEOpFromTFEOp)
  get_diagopz(odeop.tfeop)
end

function get_Tfunc(odeop::TSWODEOpFromTFEOp)
  get_Tfunc(odeop.tfeop)
end

function get_Tspace(odeop::TSWODEOpFromTFEOp)
  get_Tspace(odeop.tfeop)
end


function get_bbarfunc(odeop::TSWODEOpFromTFEOp)
  get_bbarfunc(odeop.tfeop)
end

function get_bbarspace(odeop::TSWODEOpFromTFEOp)
  get_bbarspace(odeop.tfeop)
end


function is_const_jac(odeop::TSWODEOpFromTFEOp)
  is_const_jac(odeop.tfeop)
end


function Gridap.ODEs.allocate_odeopcache(
  odeop::TSWODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  ys,b,w,z
)
  # println("allocating TSW ode op cache")

  progop = get_progop(odeop) # u, h, B
  diagop = get_diagop(odeop) # F, Φ, q
  diagopb = get_diagopb(odeop) # b
  diagopw = get_diagopw(odeop) # bbar
  diagopz = get_diagopz(odeop) # btilde

  ts = (0,0.01)

  # Allocate FE spaces for derivatives
  order = Gridap.ODEs.get_order(odeop)
  Ut = Gridap.ODEs.get_trial(progop)
  U = Gridap.ODEs.allocate_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for k in 1:order
    Uts = (Uts..., ∂t(Uts[k]))
    Us = (Us..., Gridap.ODEs.allocate_space(Uts[k+1]))
  end


  # Allocate the cache of the FE operator
  tfeopcache = allocate_tfeopcache(odeop.tfeop, t, us) # nothing

  # Variables for assembly
  uh = FEFunction(Us[1], us[1]) #_make_uh_from_us(odeop, us, Us)





  ############# DIAGNOSTICS ####################

  # Allocate FE spaces for derivatives - diagnostics F, Φ, q
  order = get_order(odeop)
  Yt = get_trial(diagop)
  Y = allocate_space(Yt)
  Yts = (Yt,)
  Ys = (Y,)
  for k in 1:order
    Yts = (Yts..., ∂t(Yts[k]))
    Ys = (Ys..., allocate_space(Yts[k+1]))
  end

  # set initial diagnostics to be zero
  yh = FEFunction(Ys[1],ys)
  diagslvrcache = nothing

  # Allocate FE spaces for derivatives - diagnostics b
  Bt = get_trial(diagopb)
  B = allocate_space(Bt)
  Bts = (Bt,)
  Bs = (B,)
  for k in 1:order
    Bts = (Bts..., ∂t(Bts[k]))
    Bs = (Bs..., allocate_space(Bts[k+1]))
  end

  # set initial diagnostics to be zero
  bh = FEFunction(Bs[1],b)

  # allocate in domain for diagnostics
  _Yt = evaluate(Gridap.FESpaces.get_trial(diagop), nothing)
  dy = get_trial_fe_basis(_Yt)
  V = Gridap.FESpaces.get_test(diagop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagop)

  jac_diag = Gridap.ODEs.get_jacs(diagop)
  dc = DomainContribution()
  dc = dc + jac_diag(ts, uh, uh, yh, dy, v)
  matdata = collect_cell_matrix(_Yt, V, dc)
  J_diag = allocate_matrix(assembler, matdata)
  diagnostics = Gridap.Algebra.allocate_in_domain(J_diag)
  fill!(diagnostics,0.0)

  # println("made TSW cache")



  ## bn diagnostics
  _Bt = evaluate(Gridap.FESpaces.get_trial(diagopb), nothing)
  db = get_trial_fe_basis(_Bt)
  V = Gridap.FESpaces.get_test(diagopb)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopb)

  jac_diagb = Gridap.ODEs.get_jacs(diagopb)
  dc = DomainContribution()
  dc = dc + jac_diagb(ts, uh, uh, bh, db, v, bh)
  matdata = collect_cell_matrix(_Bt, V, dc)
  J_diagb = allocate_matrix(assembler, matdata)
  diagnosticsb = Gridap.Algebra.allocate_in_domain(J_diagb)
  fill!(diagnosticsb,0.0)

  diagslvrcacheb = nothing

  ## bhat diagnostics -- W

  Wt = get_trial(diagopw)
  W = allocate_space(Wt)
  Wts = (Wt,)
  Ws = (W,)
  for k in 1:order
    Wts = (Wts..., ∂t(Wts[k]))
    Ws = (Ws..., allocate_space(Wts[k+1]))
  end

  # set initial diagnostics to be zero
  wh = FEFunction(Ws[1],w)

  _Wt = evaluate(Gridap.FESpaces.get_trial(diagopw), nothing)
  dw = get_trial_fe_basis(_Wt)
  V = Gridap.FESpaces.get_test(diagopw)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopw)

  jac_diagw = Gridap.ODEs.get_jacs(diagopw)
  dc = DomainContribution()
  dc = dc + jac_diagw(ts, uh, uh, wh, dw, v)
  matdata = collect_cell_matrix(_Wt, V, dc)
  J_diagw = allocate_matrix(assembler, matdata)
  diagnosticsw = Gridap.Algebra.allocate_in_domain(J_diagw)
  fill!(diagnosticsw,0.0)
  diagslvrcachew = nothing


  ## btilde diagnostics -- Z

  Zt = get_trial(diagopz)
  Z = allocate_space(Zt)
  Zts = (Zt,)
  Zs = (Z,)
  for k in 1:order
    Zts = (Zts..., ∂t(Zts[k]))
    Zs = (Zs..., allocate_space(Zts[k+1]))
  end

  # set initial diagnostics to be zero
  zh = FEFunction(Zs[1],z)

  _Zt = evaluate(Gridap.FESpaces.get_trial(diagopz), nothing)
  dz = get_trial_fe_basis(_Zt)
  V = Gridap.FESpaces.get_test(diagopz)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopz)

  bbarfunc = get_bbarfunc(odeop)
  bbarh = bbarfunc(bh,bh)

  jac_diagz = Gridap.ODEs.get_jacs(diagopz)
  dc = DomainContribution()
  dc = dc + jac_diagz(ts, yh, bbarh, zh, dz, v, bh)
  matdata = collect_cell_matrix(_Zt, V, dc)
  J_diagz = allocate_matrix(assembler, matdata)
  diagnosticsz = Gridap.Algebra.allocate_in_domain(J_diagz)
  fill!(diagnosticsz,0.0)
  diagslvrcachez = nothing

  ## Tspace

  Tt = get_Tspace(odeop)
  T = allocate_space(Tt)
  Tts = (Tt,)
  Ts = (T,)
  for k in 1:order
    Tts = (Tts..., ∂t(Tts[k]))
    Ts = (Ts..., allocate_space(Tts[k+1]))
  end

  ## bbarspace

  Bbart = get_bbarspace(odeop)
  Bbar = allocate_space(Bbart)
  Bbarts = (Bbart,)
  Bbars = (Bbar,)
  for k in 1:order
    Bbarts = (Bbarts..., ∂t(Bbarts[k]))
    Bbars = (Bbars..., allocate_space(Bbarts[k+1]))
  end


  TSWODEOpFromTFEOpCache(Us, Uts, tfeopcache,
                        Ys, Yts, diagnostics,  diagslvrcache,
                        Bs, Bts, diagnosticsb, diagslvrcacheb,
                        Ws, Wts, diagnosticsw, diagslvrcachew,
                        Zs, Zts, diagnosticsz, diagslvrcachez,
                        Ts, Tts, Bbars, Bbarts
                        )


end

function Gridap.ODEs.update_odeopcache!(odeopcache, odeop::TSWODEOpFromTFEOp, t::Real)

  Us = ()
  order = Gridap.ODEs.get_order(odeop)

  for k in 0:order
    Us = (Us..., evaluate!(odeopcache.Us[k+1], odeopcache.Uts[k+1], t))
  end
  odeopcache.Us = Us

  tfeopcache, tfeop = odeopcache.tfeopcache, odeop.tfeop
  odeopcache.tfeopcache = Gridap.ODEs.update_tfeopcache!(tfeopcache, tfeop, t)


  # update diagnostics
  Ys = ()
  for k in 0:get_order(odeop)
    Ys = (Ys..., evaluate!(odeopcache.Ys[k+1], odeopcache.Yts[k+1], t))
  end
  odeopcache.Ys = Ys

  # update diagnostics 0
  Bs = ()
  for k in 0:get_order(odeop)
    Bs = (Bs..., evaluate!(odeopcache.Bs[k+1], odeopcache.Bts[k+1], t))
  end
  odeopcache.Bs = Bs

  # update diagnostics W
  Ws = ()
  for k in 0:get_order(odeop)
    Ws = (Ws..., evaluate!(odeopcache.Ws[k+1], odeopcache.Wts[k+1], t))
  end
  odeopcache.Ws = Ws

  # update diagnostics Z
  Zs = ()
  for k in 0:get_order(odeop)
    Zs = (Zs..., evaluate!(odeopcache.Zs[k+1], odeopcache.Zts[k+1], t))
  end
  odeopcache.Zs = Zs

  # update T
  Ts = ()
  for k in 0:get_order(odeop)
    Ts = (Ts..., evaluate!(odeopcache.Ts[k+1], odeopcache.Tts[k+1], t))
  end
  odeopcache.Ts = Ts

  # update bbar
  Bbars = ()
  for k in 0:get_order(odeop)
    Bbars = (Bbars..., evaluate!(odeopcache.Bbars[k+1], odeopcache.Bbarts[k+1], t))
  end
  odeopcache.Bbars = Bbars


  odeopcache


end


function copy_to_cache!(odeopcache, b,y,w,z)
  copy!(odeopcache.diagnosticsb, b)
  copy!(odeopcache.diagnostics, y)
  copy!(odeopcache.diagnosticsw, w)
  copy!(odeopcache.diagnosticsz, z)

  odeopcache

end



function update_diagnostics!(odeopcache, y, cache)

  # copy!(odeopcache.diagnostics, y)
  odeopcache.diagslvrcache = cache
  odeopcache

end

function update_diagnosticsb!(odeopcache, b, cache)

  # copy!(odeopcache.diagnostics0, y0)
  odeopcache.diagslvrcacheb = cache
  odeopcache

end


function update_diagnosticsw!(odeopcache, w, cache)

  # copy!(odeopcache.diagnosticsw, w)
  odeopcache.diagslvrcachew = cache
  odeopcache

end


function update_diagnosticsz!(odeopcache, z, cache)

  # copy!(odeopcache.diagnosticsz, z)
  odeopcache.diagslvrcachez = cache
  odeopcache

end

function allocate_prognostic_jacobian(
  odeop::TSWODEOpFromTFEOp,
  ts::Tuple, u,
  odeopcache, u0,
  b, y, w, z, b0
)
  progop = get_progop(odeop)

  Us = odeopcache.Us
  Bs = odeopcache.Bs
  Ys = odeopcache.Ys
  Zs = odeopcache.Zs
  Bbars = odeopcache.Bbars

  uh0 = EvaluationFunction(Us[1], u0)
  uh = EvaluationFunction(Us[1], u)
  bh = EvaluationFunction(Bs[1],b)
  yh =  EvaluationFunction(Ys[1],y)
  zh = EvaluationFunction(Zs[1],z)
  bh0 = FEFunction(Bs[1], b0)
  bbarfunc = get_bbarfunc(odeop)
  bbarh = interpolate(bbarfunc(bh0,bh), Bbars[1])

  V = Gridap.ODEs.get_test(progop)
  v = get_fe_basis(V)
  Ut = evaluate(Gridap.ODEs.get_trial(progop), nothing)
  du = get_trial_fe_basis(Ut)
  assembler = Gridap.ODEs.get_assembler(progop)

  #### Compute and store Prognostic jacobian
  jac_prog = Gridap.ODEs.get_jacs(progop)
  dc = DomainContribution()
  dc = dc +  jac_prog(ts, uh0, uh, du, v,bh,yh, zh, bbarh)
  matdata = collect_cell_matrix(Ut, V, dc)
  allocate_matrix(assembler, matdata)

end


function prognostic_jacobian!(J, odeop::TSWODEOpFromTFEOp, ts::Tuple, u, odeopcache,
        u0, b, y, w, z,b0)
  # println("updating prognostic jacobian")

  ## update prognostic jacobian
  progop = get_progop(odeop)

  Us = odeopcache.Us
  Bs = odeopcache.Bs
  Ys = odeopcache.Ys
  Zs = odeopcache.Zs
  Bbars = odeopcache.Bbars

  uh0 = EvaluationFunction(Us[1], u0)
  uh = EvaluationFunction(Us[1], u)
  bh = EvaluationFunction(Bs[1], b)
  yh = EvaluationFunction(Ys[1], y)
  zh = EvaluationFunction(Zs[1], z)
  bh0 = FEFunction(Bs[1], b0)
  bbarfunc = get_bbarfunc(odeop)
  bbarh = interpolate(bbarfunc(bh0,bh), Bbars[1])

  Ut = evaluate(Gridap.FESpaces.get_trial(progop), nothing)
  du = get_trial_fe_basis(Ut)
  V = Gridap.FESpaces.get_test(progop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(progop)

  jac = Gridap.ODEs.get_jacs(progop)
  dc = DomainContribution()
  dc = dc + jac(ts, uh0, uh, du, v, bh, yh, zh, bbarh)

  matdata = collect_cell_matrix(Ut, V, dc)
  assemble_matrix!(J, assembler, matdata)

  J

end
