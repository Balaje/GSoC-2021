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
model = simplexify(CartesianDiscreteModel(domain,partition))
btrian = BoundaryTriangulation(model,tags=[7,8])

GridTopology(btrian)
