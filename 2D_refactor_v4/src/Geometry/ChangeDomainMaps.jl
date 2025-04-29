"""
Let Ωp be the parametric space of central angles α,β on panel 1
Let Ωa be the ambient space oon the sphere X_s, Y_s, Z_s

Let M: Ωp → Ωa be the map from the parametric space on panel 1 to the ambient
space on the sphere, M = R_1p ∘ σ ∘ γ, so that
  (X_s,Y_s,Z_s) = M(α,β)
                = R_1p( σ( γ( α,β) ) )

Then Minv: Ωa → Ωp is the map from the ambient space on the sphere to the
parametric space on panel 1, Minv = γinv ∘ σinv ∘ R_p1, so that
  (α,β) = Minv(X_s,Y_s,Z_s)
        = γinv ( σinv ( R_p1(X_s,Y_s,Z_s) ) )

Consider a function f: Ωp → R. Then Minv ∘ f: (Ωa → Ωp) ∘ (Ωp → R) = Ωa → R
Consider a function g: Ωa → R. Then    M ∘ f: (Ωp → Ωa) ∘ (Ωa → R) = Ωp → R

"""
# R1pFields = lazy_map(x-> PanelRotationField(r1p[x]), 1:6)
# Rp1Fields = lazy_map(x-> PanelRotationField(rp1[x]), 1:6)
γ = GnomonicField()
γinv = InvGnomonicField()
σ = SigmaField(r) # sigma is the inverse of itself

M = lazy_map(x->    PanelRotationField(r1p[x]) ∘ σ ∘ γ, 1:6)      # parametric -> ambient
Minv = lazy_map(x-> γinv ∘ σ ∘ PanelRotationField(rp1[x]), 1:6)  # ambient -> parametric

L = lazy_map(x->   σ ∘ PanelRotationField(r1p[x]) ∘ σ ∘ γ, 1:6)      # parametric -> ambient -> latlon_p
Linv = lazy_map(x-> γinv ∘ σ ∘ PanelRotationField(rp1[x]) ∘ σ, 1:6) # latlon_p -> ambient -> parametric
