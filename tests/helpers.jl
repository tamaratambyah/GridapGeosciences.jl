l2(e,dΩ) = sum(∫( e⋅e )dΩ)
nc(panel_model) = num_cells(panel_model)/6 ## nc = num cells per panel
dx(nc) = sqrt( 4*π*RADIUS^2 / (6*sqrt(nc)^2) )
nref(nc) = Int(log2(sqrt(nc))) ## level of refinement
