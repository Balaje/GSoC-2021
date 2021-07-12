using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Arrays
using Gridap
using Test


# Create domain
domain = (0,4,0,4)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
labels = get_face_labeling(model)
bgface_to_mask = get_face_mask(labels,[5,6,7,8],1)
bmodel = BoundaryDiscreteModel(Polytope{1},model,bgface_to_mask)

reffe = ReferenceFE(lagrangian, Float64, 1)
f(x) = x[1] + x[2]
V1 = FESpace(bmodel, reffe)
fh = interpolate_everywhere(f, V1)

# Works if gridtopology.jl is included.
pt = Point(rand(2))
@test evaluate(fh, pt) ≈ evaluate(f, pt)



#GridTopology(btrian)

#GridTopology(bmodel)
