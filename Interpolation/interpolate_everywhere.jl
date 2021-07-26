using Gridap
using Gridap.ReferenceFEs
using Gridap.Arrays
import Gridap.Arrays.evaluate!
import Gridap.CellData.CellField
using Test

function CellField(V::SingleFieldFESpace, fh::FEFunction)
  trian = get_triangulation(V)
  b = get_fe_dof_basis(V)
  cell_vals = _some_function_to_be_named(b, fh)
  CellField(V, cell_vals)
end

function _some_function_to_be_named(b::CellDof, fh)
  bs = get_data(b)
  trian = get_triangulation(b)
  cache = return_cache(testitem(bs), fh)
  cell_coords = get_cell_points(b).cell_phys_point
  cell_vals = lazy_map(i -> evaluate!(cache, bs[i], fh, cell_coords[i]), 1:num_cells(trian))
end

function evaluate!(cache, b::MomentBasedDofBasis, field, points)
  c, cf = cache
  vals = evaluate!(cf, field, points)
  dofs = c.array
  ReferenceFEs._eval_moment_dof_basis!(dofs, vals, b)
  dofs
end

function evaluate!(cache, b::LagrangianDofBasis, field, points)
  c, cf = cache
  vals = evaluate!(cf,field,points)
  ndofs = length(b.dof_to_node)
  T = eltype(vals)
  ncomps = num_components(T)
  ReferenceFEs._evaluate_lagr_dof!(c,vals,b.node_and_comp_to_dof,ndofs,ncomps)
end

function interpolate_everywhere_non_compatible_trian(fh::FEFunction, V::SingleFieldFESpace)
  cell_field = CellField(V, fh)
  gh = interpolate_everywhere(cell_field, V)
end

p = QUAD
D = num_dims(QUAD)
et = Float64
source_model = CartesianDiscreteModel((0,1,0,1),(10,10))
pt = VectorValue(rand(2))

@testset "Test interpolation Lagrangian" begin
  # Lagrangian space -> Lagrangian space
  f(x) = x[1] + x[2]
  reffe = LagrangianRefFE(et, p, 1)
  V₁ = FESpace(source_model, reffe, conformity=:H1)
  fh = interpolate_everywhere(f, V₁)
  # Target Lagrangian Space
  reffe = LagrangianRefFE(et, p, 2)
  model = CartesianDiscreteModel((0,1,0,1),(40,40))
  V₂ = FESpace(model, reffe, conformity=:H1)
  gh = interpolate_everywhere_non_compatible_trian(fh, V₂)

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
  reffe = RaviartThomasRefFE(et, p, 2)
  model = CartesianDiscreteModel((0,1,0,1),(40,40))
  V₂ = FESpace(model, reffe, conformity=:HDiv)
  gh = interpolate_everywhere_non_compatible_trian(fh, V₂)

  pts = [VectorValue(rand(2)) for i=1:10]
  for pt in pts
    @test gh(pt) ≈ fh(pt)
  end
end
