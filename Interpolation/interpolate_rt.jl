using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap

p = QUAD
D = num_dims(QUAD)
et = Float64
f(x) = VectorValue([x[1]+x[2],0])
reffe = RaviartThomasRefFE(et, p, 1)

# Source RT Space.
model = CartesianDiscreteModel((0,1,0,1),(10,10))
V1 = FESpace(model, reffe, conformity=:HDiv)
fh = interpolate_everywhere(f, V1)

@test evaluate(fh, Point(0.1,0.1)) ≈ evaluate(f, Point(0.1,0.1))
