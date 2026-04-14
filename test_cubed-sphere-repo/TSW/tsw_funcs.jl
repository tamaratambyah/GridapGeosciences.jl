function get_FEspaces(model,p)
  # u space, F space
  V = TestFESpace(model,
                    ReferenceFE(raviart_thomas,Float64,p),
                    conformity=:Hdiv
                    )
  U = TransientTrialFESpace(V)

  # h space, B space, Φ space, T space, b space
  W = TestFESpace(model,
                  ReferenceFE(lagrangian ,Float64,p),
                  conformity=:L2)
  R = TransientTrialFESpace(W)

  # b^2 space
  W_2 = TestFESpace(model,
                  ReferenceFE(lagrangian ,Float64,2*(p)),
                  conformity=:L2)
  R_2 = TransientTrialFESpace(W_2)

  # q space
  P = TestFESpace(model,
                  ReferenceFE(lagrangian,Float64,p+1),
                  conformity=:H1)
  H = TransientTrialFESpace(P)

  spaces = (;U=U,V=V,R=R,W=W,R_2=R_2,W_2=W_2,H=H,P=P )

  X_prog = TransientMultiFieldFESpace([U,R,R]) # u, h, B
  Y_prog = MultiFieldFESpace([V,W,W]) # u, h, B

  X_diag = TransientMultiFieldFESpace([U,R,H,H]) # F, Φ, q, ω
  Y_diag = MultiFieldFESpace([V,W,P,P]) # F, Φ, q, ω

  multi_spaces = (;X_prog=X_prog,Y_prog=Y_prog,X_diag=X_diag,Y_diag=Y_diag)

  return spaces, multi_spaces

end

function compute_casimirs(u,h,B,b,b3,F,dΩ,dΛ,n_Λ,b0,upwinding,ε,soft, h0,B0,dt)

  energy =  sum( ∫(0.5*h*(u⋅(metric_cf⋅u))*meas_cf + 0.5*h*B*meas_cf  )dΩ)
  entropy = sum( ∫( 0.5*b*b*h*meas_cf )dΩ  )
  _entropy = sum( ∫( 0.5*b*b0*h*meas_cf )dΩ  )
  entropy2 = sum( ∫( 0.5*B*B/h*meas_cf )dΩ  )

  ##### TO DO:::: PULLBACKS BELOW !!!!
  entropy_internal_btilde = sum( ∫( -1.0*(∇⋅F)*(0.5)*( b0*b0 + b*b )  )dΩ
                              + ∫( (b3*(0.5*(b0+b))  )*(∇⋅F) )dΩ  )

  entropy_internal2 = sum( ∫( 0.5*(0.5*(b0 + b))*(F⋅∇((0.5*(b0 + b))))  )dΩ - ∫( 0.5*(∇((0.5*(b0 + b)))⋅F)*(0.5*(b0 + b)) )dΩ )

  entropy_center = sum( ∫( 0.5*(mean(F*(0.5*(b0 + b)))⋅jump((0.5*(b0 + b))*n_Λ ) ) )dΛ - ∫( 0.5*(mean(F*(0.5*(b0 + b)))⋅jump((0.5*(b0 + b))*n_Λ ))  )dΛ   )

  entropy_upwind = 0.0

  time_error = sum( ∫( 0.5*( b0*h0 + b*h )*(b-b0)/dt)dΩ
                   -∫( 0.5*(b0 + b)*(B-B0)/dt)dΩ
                   +∫( 0.5*( b0*b0 + b*b )*(h-h0)/dt )dΩ)

  mass = sum( ∫( h )dΩ)
  denB = sum( ∫( h*b )dΩ)

  if upwinding
    entropy_upwind = sum( ∫( (upwinding_sign∘( (F⋅n_Λ).plus, ε, soft )*( (F⋅n_Λ).plus ) )*jump((0.5*(b0 + b))*n_Λ )⋅jump((0.5*(b0 + b))*n_Λ )   )dΛ )
  end

  casimirs = (;E=energy,S=entropy,S2=entropy2,Sl=_entropy,
              Sinternal_btilde=entropy_internal_btilde,Sinternal2=entropy_internal2,
              Scenter=entropy_center,Supwind=entropy_upwind,
              time_error=time_error,
              M=mass,denB=denB)
  return casimirs
end


Tfunc((u0,h0,B0),(u,h,B)) = 0.25*(h + h0 )
bbarfunc(b0,b) = 0.5*( b0 + b )


########################################
# Diagnostics # bhat (b2)
# diagnosed after each bbar compute
########################################

function res_diagw(dΩ)
  _res_diagw((t,dt),(u0,h0,B0),(u,h,B),(b0,λ0),(b,λ),b2,l2) = (
          ∫( b2*l2*meas_cf - (1.0/2.0)*( b0*b0 + b*b )*l2*meas_cf   )dΩ
  )
end

function jac_diagw(dΩ)
  _jac_diagw((t,dt),(u0,h0,B0),(u,h,B),b2,db2,l2) = (
          ∫( db2*l2*meas_cf  )dΩ
  )
end
