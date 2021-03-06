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
using Test

domain = (0,1,0,1)
partition = (20,20)
# Old reffe (order 1)
reffe = ReferenceFE(lagrangian, Float64, 1)

# Construct a solution in the Cartesian Space.
model = CartesianDiscreteModel(domain, partition)
Ω1 = Triangulation(model)
V1 = FESpace(model, reffe)
f(x) = x[1] + x[2]
fh = interpolate_everywhere(f, V1) # Old FEFunction
writevtk(Ω1,"source",cellfields=["fh"=>fh])

# Second FE Space
function rndm(p::Point)
  r, s = p
  x = r + 0.05*sin(2π*r)*sin(2π*s)
  y = s + 0.05*sin(2π*r)*sin(2π*s)
  Point(x,y)
end
partition=(40,40)
model2 = CartesianDiscreteModel(domain,partition; map=rndm)
Ω2 = Triangulation(model2)
# New reffe (order 2)
reffe = ReferenceFE(lagrangian, Float64, 2)
V2 = FESpace(model2, reffe)

gh = interpolate_everywhere_non_compatible_trian(fh, V2)


# Test evaluate on gh (Generic Cell Field): Fails for some points
# Modify nlsolve in src/Fields/InverseFields.jl
#     res = nlsolve(df,y₀,method=:newton,linesearch=BackTracking())
# to make everything work
@show evaluate(gh, Point(0.1,0.1)) # --Works
@show evaluate(gh, Point(0.1,0.5)) # --Does not work
@show evaluate(gh, VectorValue(rand(2))) # --Does not work
#evaluate(gh, Point(0.1,0.56)) # --Does not work
# Works sometimes ...
for i ∈ 1:10
    x = VectorValue(rand(2))
    @test evaluate(gh, x) ≈ f(x)
end
