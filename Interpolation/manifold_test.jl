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
btrian = BoundaryTriangulation(model,tags=[5,6,7,8])

#GridTopology(btrian)

labels = get_face_labeling(model)
bgface_to_mask = get_face_mask(labels,[5,6,7,8],1)
bmodel = BoundaryDiscreteModel(Polytope{1},model,bgface_to_mask)

GridTopology(bmodel)
