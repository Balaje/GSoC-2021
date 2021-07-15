using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Arrays
using Gridap
using Test


# Create domain
domain = (0,1,0,1,0,1)
cells = (3,3,3)
bgmodel = simplexify(CartesianDiscreteModel(domain,cells))
#bgmodel = CartesianDiscreteModel(domain,cells)
labels = get_face_labeling(bgmodel)
bgface_to_mask = get_face_mask(labels,[21,22,23,24],2)
model = BoundaryDiscreteModel(Polytope{2},bgmodel,bgface_to_mask)

reffe = ReferenceFE(lagrangian, Float64, 2)
V₁ = FESpace(model, reffe)
@show V₁
f(x) = x[1]^2+x[2]
fh = interpolate_everywhere(f, V₁)
