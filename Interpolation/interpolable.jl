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

function _old_cell_vals(V::SingleFieldFESpace, f::Interpolatable)
  fe_basis = get_fe_dof_basis(V)
  trian = get_triangulation(V)
  fe_basis_phys = change_domain(fe_basis, ReferenceDomain(), PhysicalDomain())
  bs = get_data(fe_basis_phys)
  cache = return_cache(testitem(bs), f.uh)
  cell_vals = lazy_map(i -> evaluate!(cache, bs[i], f.uh), 1:num_cells(trian))
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

ifh = Interpolatable(fh)

# Check whether the new implementation is the same as old
print("Check whether: ")
@show FESpaces._cell_vals(V₂, ifh) == _old_cell_vals(V₂, ifh)
print("\nNew implemenation\n")
@btime FESpaces._cell_vals(V₂, ifh);
print("\nOld implemenation\n")
@btime _old_cell_vals(V₂, ifh);
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

@testset "Test interpolation RT" begin
  gh = interpolate_everywhere(ifh, V₂)
  pts = [VectorValue(rand(2)) for i=1:10]
  for pt in pts
    @test gh(pt) ≈ fh(pt)
  end
end
