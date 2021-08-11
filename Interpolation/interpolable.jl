#include("interpolable_include.jl")

using Gridap
using Test
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using Gridap.Geometry
using StaticArrays
using NearestNeighbors
using Gridap.Helpers

# struct Interpolable{A} <: Function
#   uh::A
#   tol::Float64
#   cache
#   searchmethod
#   function Interpolable(uh, cache; tol=1e-6, searchmethod=:kdtree)
#     new{typeof(uh)}(uh, tol, cache, searchmethod)
#   end
# end

# function Interpolable(uh::CellField; tol=1e-6, searchmethod=:kdtree)
#   if(searchmethod != :kdtree)
#     @notimplemented
#   end
#   trian = get_triangulation(uh)
#   cache1 = CellData._point_to_cell_cache(trian)

#   cell_f = get_array(uh)
#   cell_f_cache = array_cache(cell_f)
#   cf = testitem(cell_f)
#   T = eltype(testitem(trian.node_coords))
#   dim = num_point_dims(trian)
#   f_cache = return_cache(cf,VectorValue(rand(T,dim)))
#   cache2 = cell_f_cache, f_cache, cell_f, uh
#   cache = cache1, cache2

#   Interpolable(uh, cache; tol, searchmethod)
# end

# Arrays.return_cache(f::Interpolable, x::Point) = f.cache
# Arrays.return_cache(f::Interpolable, xs::AbstractVector{<:Point}) = f.cache
# Arrays.evaluate!(cache, f::Interpolable, x::Point) = evaluate!(cache, f.uh, x)
# Arrays.evaluate!(cache, f::Interpolable, xs::AbstractVector{<:Point}) = evaluate!(cache, f.uh, xs)

"""
 interpolate_everywhere(Interpolatable(fₕ), V₂)  Works without these routines.
   Here fₕ ∈ V₁ ≠ V₂
 Check interpolable_1.jl

Functions for Optimizing interpolations:
- _cell_vals
"""

# function FESpaces._cell_vals(V::SingleFieldFESpace, f::Interpolable)
#   fe_basis = get_fe_dof_basis(V)
#   trian = get_triangulation(V)
#   cache = return_cache(f, Point(rand(2)))
#   b = change_domain(fe_basis, ReferenceDomain(), PhysicalDomain())
#   cf = CellField(x->evaluate!(cache, f, x), trian, ReferenceDomain())
#   lazy_map(evaluate, get_data(b), get_data(cf))
# end

"""
Some Tests with optimizations...
"""

using BenchmarkTools

p = QUAD
D = num_dims(QUAD)
et = Float64
source_model = CartesianDiscreteModel((0,1,0,1),(2,2))

@testset "Test interpolation Lagrangian" begin
  f(x) = x[1] + x[2]
  reffe = LagrangianRefFE(et, p, 1)
  V₁ = FESpace(source_model, reffe, conformity=:H1)
  fh = interpolate_everywhere(f, V₁)
  # Target Lagrangian Space
  reffe = LagrangianRefFE(et, p, 2)
  model = CartesianDiscreteModel((0,1,0,1),(4,4))
  V₂ = FESpace(model, reffe, conformity=:H1)

  ifh = Interpolable(fh)
  # Lagrangian space -> Lagrangian space
  try
    interpolate_everywhere(fh, V₂)
  catch
    gh = interpolate_everywhere(ifh, V₂)
    pts = [VectorValue(rand(2)) for i=1:10]
    for pt in pts
      @test gh(pt) ≈ fh(pt)
    end
  end

  # VectorValued Lagrangian
  fᵥ(x) = VectorValue([x[1], x[1]+x[2]])
  reffe = ReferenceFE(lagrangian, VectorValue{2,et}, 1)
  V₁ = FESpace(source_model, reffe, conformity=:H1)
  fh = interpolate_everywhere(fᵥ, V₁)
  # Target
  reffe = ReferenceFE(lagrangian, VectorValue{2,et},  2)
  V₂ = FESpace(model, reffe, conformity=:H1)

  ifh = Interpolable(fh);
  gh = interpolate_everywhere(ifh, V₂)
  pts = [VectorValue(rand(2)) for i=1:10]
  for pt in pts
    @test gh(pt) ≈ fh(pt)
  end
end

@testset "Test interpolation RT" begin
  # RT Space -> RT Space
  f(x) = VectorValue([x[1], x[2]])
  reffe = RaviartThomasRefFE(et, p, 0)
  V₁ = FESpace(source_model, reffe, conformity=:HDiv)
  fh = interpolate_everywhere(f, V₁);
  # Target RT Space
  reffe = RaviartThomasRefFE(et, p, 1)
  model = CartesianDiscreteModel((0,1,0,1),(40,40))
  V₂ = FESpace(model, reffe, conformity=:HDiv)

  ifh = Interpolable(fh)
  gh = interpolate_everywhere(ifh, V₂)
  pts = [VectorValue(rand(2)) for i=1:10]
  for pt in pts
    @test gh(pt) ≈ fh(pt)
  end
end

@testset "Test interpolation Multifield" begin
  f₁(x) = x[1]+x[2]
  f₂(x) = x[1]
  # Source FESpace
  reffe = LagrangianRefFE(et, p, 1)
  V₁ = FESpace(source_model, reffe, conformity=:H1)
  V₁² = MultiFieldFESpace([V₁,V₁])
  fh = interpolate_everywhere([f₁, f₂], V₁²)

  # Target Lagrangian FESpace
  reffe = LagrangianRefFE(et, p, 2)
  model = CartesianDiscreteModel((0,1,0,1), (40,40))
  V₂ = FESpace(model, reffe, conformity=:H1)
  V₂² = MultiFieldFESpace([V₂,V₂])

  fh₁,fh₂ = fh
  ifh₁ = Interpolable(fh₁)
  ifh₂ = Interpolable(fh₂)

  gh = interpolate_everywhere([ifh₁,ifh₂], V₂²)

  pts = [VectorValue(rand(2)) for i=1:10]
  gh₁,gh₂ = gh
  for pt in pts
    @test gh₁(pt) ≈ fh₁(pt)
    @test gh₂(pt) ≈ fh₂(pt)
  end
end
