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
function rndm(p::Point)
    r, s = p
    x = r + 0.02*sin(2π*r)*sin(2π*s)
    y = s + 0.02*sin(2π*r)*sin(2π*s)
    Point(x,y)
end
partition=(40,40)
model2 = CartesianDiscreteModel(domain,partition; map=rndm)
Ω2 = Triangulation(model2)
reffe = ReferenceFE(lagrangian, Float64, 2)
V2 = FESpace(model2, reffe)

# Main Solution
phys_point = get_cell_points(get_fe_dof_basis(V2)).cell_phys_point
fh_phys_coords(x) = evaluate(fh, x)
gh = CellField( V2, lazy_map(fh_phys_coords, phys_point) )

# Write the solution
writevtk(get_triangulation(gh),"target_sin",cellfields=["fh"=>gh])
