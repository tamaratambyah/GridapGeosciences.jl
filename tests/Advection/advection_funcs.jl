#### Replicate test in Section 5.4 of Rognes2013 paper
vecX(XYZ) = VectorValue(-XYZ[2],XYZ[1],0.0)
u0(XYZ) = exp(-(XYZ[2]^2 + XYZ[3]^2)  )
u0vecX(XYZ) = u0(XYZ)*vecX(XYZ)
