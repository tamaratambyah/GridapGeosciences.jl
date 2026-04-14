
struct TSWGenericTransientFESolution <: TransientFESolution
  odesltn::ODESolution
  trial_prog
  trial_diag0
  trial_diag
  trial_diagw
  trial_diagz
end


# Constructors
function Gridap.ODEs.GenericTransientFESolution(
  odeslvr::TSWPoissonIntegrator, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uhs0::Tuple{Vararg{CellField}},
  bhs, yhs, whs, zhs
)
  # println("my func")



  odeop = Gridap.ODEs.get_algebraic_operator(tfeop)
  us0 = get_free_dof_values.(uhs0)
  bs = get_free_dof_values.(bhs)
  ys = get_free_dof_values.(yhs)
  ws = get_free_dof_values.(whs)
  zs = get_free_dof_values.(zhs)
  odesltn = solve(odeslvr, odeop, t0, tF, us0,bs,ys,ws,zs)
  trial_prog = Gridap.FESpaces.get_trial(tfeop )
  # trial_diag = get_diag_trial(tfeop)
  TSWGenericTransientFESolution(odesltn, trial_prog,trial_prog,trial_prog,trial_prog,trial_prog)
end

function Gridap.ODEs.GenericTransientFESolution(
  odeslvr::TSWPoissonIntegrator, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uh0::CellField,
  bh, yh, wh, zh
)
  uhs0 = (uh0,)
  bhs = (bh,)
  yhs = (yh, )
  whs = (wh, )
  zhs = (zh, )
  GenericTransientFESolution(odeslvr, tfeop, t0, tF, uhs0, bhs, yhs, whs, zhs)
end

function Base.iterate(tfesltn::TSWGenericTransientFESolution)
  ode_it = iterate(tfesltn.odesltn)
  if isnothing(ode_it)
    return nothing
  end

  ode_it_data, ode_it_state = ode_it
  tF, uF = ode_it_data

  Uh = Gridap.ODEs.allocate_space(tfesltn.trial_prog)
  Uh = evaluate!(Uh, tfesltn.trial_prog, tF)
  uhF = FEFunction(Uh, uF)



  tfe_it_data = (tF, uhF)

  tfe_it_state = (Uh, ode_it_state)
  (tfe_it_data, tfe_it_state)
end

function Base.iterate(tfesltn::TSWGenericTransientFESolution, state)
  Uh, ode_it_state = state

  ode_it = iterate(tfesltn.odesltn, ode_it_state)
  if isnothing(ode_it)
    return nothing
  end

  ode_it_data, ode_it_state = ode_it
  tF, uF = ode_it_data

  Uh = evaluate!(Uh, tfesltn.trial_prog, tF)
  uhF = FEFunction(Uh, uF)

  tfe_it_data = (tF, uhF)
  tfe_it_state = (Uh, ode_it_state)
  (tfe_it_data, tfe_it_state)
end


function Gridap.Algebra.solve(
  odeslvr::TSWPoissonIntegrator, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uhs0::Tuple{Vararg{CellField}},
  bhs, yhs, whs, zhs
)
  GenericTransientFESolution(odeslvr, tfeop, t0, tF, uhs0, bhs, yhs, whs, zhs)
end

function Gridap.Algebra.solve(
  odeslvr::TSWPoissonIntegrator, tfeop::TransientFEOperator,
  t0::Real, tF::Real, uh0::CellField,
  bh, yh, wh, zh
)
  uhs0 = (uh0,)
  bhs = (bh,)
  yhs = (yh, )
  whs = (wh, )
  zhs = (zh, )
  solve(odeslvr, tfeop, t0, tF, uhs0, bhs, yhs, whs, zhs)
end
