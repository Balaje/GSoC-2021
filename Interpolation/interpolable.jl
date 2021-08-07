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

include("interpolable_include.jl")

struct Interpolatable{A} <: Function
  uh::A
  tol::Float64
  cache
  function Interpolatable(uh; tol=1e-6)
    trian = get_triangulation(uh)
    topo = GridTopology(trian)
    vertex_coordinates = Geometry.get_vertex_coordinates(topo)
    kdtree = KDTree(map(nc -> SVector(Tuple(nc)), vertex_coordinates))
    D = num_cell_dims(trian)
    vertex_to_cells = get_faces(topo, 0, D)
    cell_to_ctype = get_cell_type(trian)
    ctype_to_reffe = get_reffes(trian)
    ctype_to_polytope = map(get_polytope, ctype_to_reffe)
    cell_map = get_cell_map(trian)
    cache1 = kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map
    new{typeof(uh)}(uh, tol, cache1)
  end
end

function Arrays.return_cache(f::Interpolatable, x::Point)
  cache1 = f.cache
  cell_f = get_array(f.uh)
  cell_f_cache = array_cache(cell_f)
  cf = testitem(cell_f)
  f_cache = return_cache(cf,x)
  cache2 = cell_f_cache, f_cache, cell_f, f.uh
  cache1, cache2
end
Arrays.return_cache(f::Interpolatable, xs::AbstractVector{<:Point}) = return_cache(f, testitem(xs))
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
  fe_basis_phys = change_domain(fe_basis, ReferenceDomain(), PhysicalDomain())
  cache = return_cache(f, Point(0,0));
  cf = CellField(x->evaluate!(cache, f, x), trian, PhysicalDomain())
  fe_basis(cf)
end


"""
Some Tests with optimizations...
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

# ifh = Interpolatable(fh)
# @testset "Test interpolation Lagrangian" begin
#   # Lagrangian space -> Lagrangian space
#   try
#     interpolate_everywhere(fh, V₂)
#   catch
#     @btime interpolate_everywhere(ifh, V₂) # Check time for interpolation
#     gh = interpolate_everywhere(ifh, V₂)
#     pts = [VectorValue(rand(2)) for i=1:10]
#     for pt in pts
#       @test gh(pt) ≈ fh(pt)
#     end
#   end

# end
