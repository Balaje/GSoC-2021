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

struct Interpolatable{A} <: Function
  uh::A
  tol::Float64
  cache
  searchmethod
  function Interpolatable(uh, cache; tol=1e-6, searchmethod=:kdtree)
    new{typeof(uh)}(uh, tol, cache, searchmethod)
  end
end

function Interpolatable(uh::CellField; tol=1e-6, searchmethod=:kdtree)
  if(searchmethod != :kdtree)
    @notimplemented
  end
  trian = get_triangulation(uh)
  cache1 = CellData._point_to_cell_cache(trian)

  cell_f = get_array(uh)
  cell_f_cache = array_cache(cell_f)
  cf = testitem(cell_f)
  T = eltype(testitem(trian.node_coords))
  dim = num_point_dims(trian)
  f_cache = return_cache(cf,VectorValue(rand(T,dim)))
  cache2 = cell_f_cache, f_cache, cell_f, uh
  cache = cache1, cache2

  Interpolatable(uh, cache; tol, searchmethod)
end

Arrays.return_cache(f::Interpolatable, x::Point) = f.cache
Arrays.return_cache(f::Interpolatable, xs::AbstractVector{<:Point}) = f.cache
Arrays.evaluate!(cache, f::Interpolatable, x::Point) = evaluate!(cache, f.uh, x)
Arrays.evaluate!(cache, f::Interpolatable, xs::AbstractVector{<:Point}) = evaluate!(cache, f.uh, xs)

"""
 interpolate_everywhere(Interpolatable(fₕ), V₂)  Works without these routines.
   Here fₕ ∈ V₁ ≠ V₂
 Check interpolable_1.jl

Functions for Optimizing interpolations:
- _cell_vals
"""

function FESpaces._cell_vals(V::SingleFieldFESpace, f::Interpolatable)
  fe_basis = get_fe_dof_basis(V)
  trian = get_triangulation(V)
  cache = return_cache(f, Point(rand(2)))
  cf = CellField(x->evaluate!(cache, f, x), trian, ReferenceDomain())
  fe_basis(cf)
end

"""
Some Tests with optimizations...
"""

using BenchmarkTools

p = QUAD
D = num_dims(QUAD)
et = Float64
source_model = CartesianDiscreteModel((0,1,0,1),(100,100))
f(x) = x[1] + x[2]
reffe = LagrangianRefFE(et, p, 1)
V₁ = FESpace(source_model, reffe, conformity=:H1)
fh = interpolate_everywhere(f, V₁)
# Target Lagrangian Space
reffe = LagrangianRefFE(et, p, 2)
model = CartesianDiscreteModel((0,1,0,1),(40,40))
V₂ = FESpace(model, reffe, conformity=:H1)

ifh = Interpolatable(fh)

print("\nTime and allocations:\n")
@btime FESpaces._cell_vals(V₂, ifh);
print("\n")

@testset "Test interpolation Lagrangian" begin
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
end

@testset "Test interpolation RT" begin
  # RT Space -> RT Space
  f(x) = VectorValue([x[1], x[2]])
  reffe = RaviartThomasRefFE(et, p, 0)
  V₁ = FESpace(source_model, reffe, conformity=:HDiv)
  fh = interpolate_everywhere(f, V₁);
  # Target RT Space
  reffe = RaviartThomasRefFE(et, p, 2)
  model = CartesianDiscreteModel((0,1,0,1),(40,40))
  V₂ = FESpace(model, reffe, conformity=:HDiv)

  ifh = Interpolatable(fh)
  gh = interpolate_everywhere(ifh, V₂)
  pts = [VectorValue(rand(2)) for i=1:10]
  for pt in pts
    @test gh(pt) ≈ fh(pt)
  end
end
