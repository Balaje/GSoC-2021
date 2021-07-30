using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Arrays
using Gridap
using Test

p = QUAD
D = num_dims(QUAD)
et = Float64
f(x) = VectorValue([x[1], x[2]])
reffe = RaviartThomasRefFE(et, p, 0)

# Source RT Space.
model = CartesianDiscreteModel((0,1,0,1),(10,10))
V1 = FESpace(model, reffe, conformity=:HDiv)
fh = interpolate_everywhere(f, V1)

# Target RT Space
reffe = RaviartThomasRefFE(et, p, 1)
model = CartesianDiscreteModel((0,1,0,1),(40,40))
V2 = FESpace(model, reffe, conformity=:HDiv)

# Now to test it
gh = interpolate_everywhere_non_compatible_trian(fh, V2)

pt = Point(rand(2))
@test evaluate(fh, pt) ≈ evaluate(gh, pt)
