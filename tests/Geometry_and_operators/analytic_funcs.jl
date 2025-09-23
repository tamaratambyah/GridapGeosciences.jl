include("../Laplace/analytic_funcs.jl")
include("../Geophysical/Williamson_functions_v2.jl")

### ambient vectors
vecX_1(XYZ) = VectorValue(-1.0*XYZ[2],XYZ[1],0.0)
vecX_2(XYZ) = VectorValue(XYZ[1]*XYZ[2],XYZ[2]*XYZ[3],XYZ[3]^2-RADIUS^2)
vecX_3(XYZ) = VectorValue(XYZ[2],XYZ[3],0.0)

ambient_vecs = Dict{Symbol,Any}()
ambient_vecs[:v1] = vecX_1
ambient_vecs[:v2] = vecX_2
ambient_vecs[:v3] = vecX_3

### williamson2 vector field - defined in latlon
williamson_vec = Dict{Symbol,Any}()
williamson_vec[:z1] = u₀(0.0)
williamson_vec[:z2] = u₀(0.05)
williamson_vec[:z3] = u₀(π/2-0.05)
williamson_vec[:z4] = u₀(π/2)
