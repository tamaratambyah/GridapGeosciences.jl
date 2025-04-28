"""
Let Ω1 be the parametric space of central angles α,β on panel 1
Let Ω2 be the ambient space oon the sphere X_s, Y_s, Z_s

Let M: Ω1 → Ω2 be the map from the parametric space on panel 1 to the ambient
space on the sphere, M = R_1p ∘ σ ∘ γ, so that
  (X_s,Y_s,Z_s) = M(α,β)
                = R_1p( σ( γ( α,β) ) )

Then Minv: Ω2 → Ω1 is the map from the ambient space on the sphere to the
parametric space on panel 1, Minv = γinv ∘ σinv ∘ R_p1, so that
  (α,β) = Minv(X_s,Y_s,Z_s)
        = γinv ( σinv ( R_p1(X_s,Y_s,Z_s) ) )

Consider a function f: Ω1 → R. Then Minv ∘ f: (Ω2 → Ω1) ∘ (Ω1 → R) = Ω2 → R
Consider a function g: Ω2 → R. Then    M ∘ f: (Ω1 → Ω2) ∘ (Ω2 → R) = Ω1 → R

"""
R1pFields = lazy_map(x-> PanelRotationField(r1p[x]), 1:6)
Rp1Fields = lazy_map(x-> PanelRotationField(rp1[x]), 1:6)
γ = GnomonicField()
γinv = InvGnomonicField()
σ = SigmaField(r) # sigma is the inverse of itself

M = lazy_map(x->    R1pFields[x]∘ σ ∘ γ, 1:6)      # parametric -> ambient
Minv = lazy_map(x-> γinv ∘ σ ∘ Rp1Fields[x], 1:6)  # ambient -> parametric
