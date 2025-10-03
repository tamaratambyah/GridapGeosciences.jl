
const options_gmres = """
-cg_ksp_type cg
-cg_ksp_rtol 1.0e-6
-cg_ksp_converged_reason
-cg_ksp_monitor
-gm_ksp_type gmres
-gm_ksp_rtol 1.0e-6
-gm_ksp_converged_reason
-gm_ksp_monitor
"""

# linear solver from options: prefix c
function petsc_ls_from_options_cg(ksp)
  @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"cg_")
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
end

# linear solver from options: prefix gm == gmres
function petsc_ls_from_options_gm(ksp)
  @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"gm_")
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
end

# linear solver - gmres, precondiioned with jacobi
function petsc_gmres_jacobi(ksp)
  rtol = PetscScalar(1.e-6)
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = GridapPETSc.PETSC.PETSC_DEFAULT

  # GMRES solver
  @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"gm_")
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)

  # Preconditioner
  pc       = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCJACOBI)

end

# linear solver - cg, precondiioned with jacobi
function petsc_cg_jacobi(ksp)
  rtol = PetscScalar(1.e-6)
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = GridapPETSc.PETSC.PETSC_DEFAULT

  # GMRES solver
  @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"cg_")
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)

  # Preconditioner
  pc       = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCJACOBI)

end
