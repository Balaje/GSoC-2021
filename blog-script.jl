using Gridap
using Gridap.Algebra
using Gridap.CellData
using Gridap.Fields
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.MultiField
using Gridap.ReferenceFEs
using StaticArrays
using Random
using Test

# Create computational domain
domain = (0,1,0,1)
partition = (5,5)
model = CartesianDiscreteModel(domain, partition)
# Sample points
pt = Point(rand(2))
pts = [Point(rand(2)) for i in 1:3]

"""
Example 1a: Evaluating FEFunction on arbitrary points. PR#523
"""
f₁(x) = x[1] + x[2]
reffe₁ = ReferenceFE(lagrangian, Float64, 1) # Lagrangian ReferenceFE
V₁ = FESpace(model, reffe₁) # Build the FESpace
fₕ = interpolate_everywhere(f₁,V₁) # Finite dimensional version of f on V₁
# Evaluating fₕ on a single point
print("\nTest evaluate Lagrange FESpaces:\n")
@show fₕ(pt)
@show f₁(pt)
@show fₕ.(pts)
@show f₁.(pts)
@show fₕ(pt) ≈ f₁(pt) # Same as evaluate(fₕ, pt) ≈ evaluate(f, pt)
@show fₕ.(pts) ≈ f₁.(pts) # Array of Points


"""
 Example 1b: Evaluating Raviart-Thomas FEFunction on arbitrary points
"""
reffe₂ = ReferenceFE(raviart_thomas, Float64, 1) # Raviart Thomas ReferenceFE
f₂(x) = VectorValue([x[1], x[2]]) # ∵ RT is a Vectorial FESpace
V₂ = FESpace(model, reffe₂)
fₕ = interpolate_everywhere(f₂, V₂) # Finite dimensional version of f on V₂
# Evaluating RT Function on arbitrary points. Fixed after PR#628
print("\nTest RT FESpaces:\n")
@show fₕ(pt)
@show f₂(pt)
@show fₕ.(pts)
@show f₂.(pts)
@show fₕ(pt) ≈ f₂(pt)
@show fₕ.(pts) ≈ f₂.(pts)


"""
Example 2a: Interpolation between Lagrange FESpaces
"""
# Define a target computational domain
domain = (0,1,0,1)
partition = (20,20)
model = CartesianDiscreteModel(domain,partition)
W₁ = FESpace(model, reffe₁)
fₕ = interpolate_everywhere(f₁,V₁)
ifₕ = Interpolable(fₕ)
gₕ = interpolate_everywhere(ifₕ,W₁) # Interpolate fₕ ∈ V₁ on to W₁
print("\nTest Interpolation on Lagrange Elements:\n")
@show gₕ(pt)
@show f₁(pt)
@show gₕ(pt) ≈ fₕ(pt) ≈ f₁(pt)


"""
Example 2b: Interpolation between RT FESpaces
"""
domain = (0,1,0,1)
partition = (20,20)
model = CartesianDiscreteModel(domain,partition)
W₂ = FESpace(model, reffe₂)
fₕ = interpolate_everywhere(f₂,V₂)
ifₕ = Interpolable(fₕ)
gₕ = interpolate_everywhere(ifₕ,W₂) # Interpolate fₕ ∈ V₂ on to W₂
print("\nTest Interpolation on RT Elements:\n")
@show gₕ(pt)
@show f₂(pt)
@show gₕ(pt) ≈ fₕ(pt) ≈ f₂(pt)
