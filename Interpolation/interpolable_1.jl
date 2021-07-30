using Gridap
using Test
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays

struct Interpolatable{A} <: Function
  uh::A
  tol::Float64
  function Interpolatable(uh; tol=1e-6)
    new{typeof(uh)}(uh,tol)
  end
end

(f::Interpolatable)(x) = f.uh(x)

"""
Some tests
"""

using BenchmarkTools

p = QUAD
D = num_dims(QUAD)
et = Float64
source_model = CartesianDiscreteModel((0,1,0,1),(10,10))
f(x) = x[1] + x[2]
reffe = LagrangianRefFE(et, p, 1)
V₁ = FESpace(source_model, reffe, conformity=:H1)
fh = interpolate_everywhere(f, V₁)
# Target Lagrangian Space
reffe = LagrangianRefFE(et, p, 2)
model = CartesianDiscreteModel((0,1,0,1),(40,40))
V₂ = FESpace(model, reffe, conformity=:H1)

# Build the Interpolatable function
ifh = Interpolatable(fh)

@testset "Test interpolation Lagrangian" begin
  # Lagrangian space -> Lagrangian space
  try
    interpolate_everywhere(fh, V₂)
  catch
    @btime interpolate_everywhere(ifh, V₂) # Check time for interpolation
    gh = interpolate_everywhere(ifh, V₂)
    pts = [VectorValue(rand(2)) for i=1:10]
    for pt in pts
      @test gh(pt) ≈ fh(pt)
    end
  end
end
