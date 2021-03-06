using Gridap
using Gridap.Algebra
using Gridap.CellData
using Gridap.Fields
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.MultiField
using Gridap.ReferenceFEs
using Gridap.Visualization
using StaticArrays
using Random

domain = (0,1,0,1)
partition = (20,20)
reffe = ReferenceFE(lagrangian, Float64, 2)

# Construct a solution in the Cartesian Space.
model = CartesianDiscreteModel(domain, partition)
Ω1 = Triangulation(model)
V1 = FESpace(model, reffe)
f(x) = x[1]^2 + x[2]^2
fh = interpolate_everywhere(f, V1) # Old FEFunction
writevtk(Ω1,"source",cellfields=["fh"=>fh])


# Second FE Space
function rndm(p::Point, D=pt_to_rand)
    r, s = p
    Random.seed!(get(D,p,1234))
    x = (r ≈ 1.0 || r ≈ 0.0) ? r : r + 0.012*rand()
    y = (s ≈ 1.0 || s ≈ 0.0) ? s : s + 0.012*rand()
    Point(x,y)
end
partition=(20,20)
# Construct a uniform partition first to build the perturbation
model = CartesianDiscreteModel(domain,partition)
# Random permutation of indices
pt_to_ind = randperm(num_nodes(model))
# Build a dictionary with random indices for points
pt_to_rand = Dict(zip(get_node_coordinates(model),pt_to_ind));
# Use the dictionary to seed points
model2 = CartesianDiscreteModel(domain,partition; map=rndm)
Ω2 = Triangulation(model2)
reffe = ReferenceFE(lagrangian, Float64, 2)
V2 = FESpace(model2, reffe)

# Main Solution
gh = interpolate_everywhere_non_compatible_trian(fh, V2)

# Test evaluate on gh (Generic Cell Field): Fails for some points
@show evaluate(gh, VectorValue(rand(2))) # Works sometimes ...
