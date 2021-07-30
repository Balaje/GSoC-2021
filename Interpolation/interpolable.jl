using Gridap
using Test
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays

struct Interpolatable <: Function
  uh::FEFunction
  tol::Float64
  function Interpolatable(uh; tol=1e-6)
    new(uh,tol)
  end
end

(f::Interpolatable)(x) = f.uh(x)

"""
 interpolate_everywhere(Interpolatable(fₕ), V₂)  Works without these routines.
   Here fₕ ∈ V₁ ≠ V₂
 Check interpolable_1.jl

Functions for Optimizing interpolations:
- _cell_vals
- _to_physical_domain
"""
function FESpaces._cell_vals(V::SingleFieldFESpace, f::Interpolatable)
  fe_basis = get_fe_dof_basis(V)
  fh = f.uh
  trian = get_triangulation(V)
  cell_maps = get_cell_map(trian)
  bs = get_data(fe_basis)
  cache = return_cache(testitem(bs), fh)
  if(DomainStyle(fe_basis) == ReferenceDomain())
    cell_vals = lazy_map(i -> evaluate!(cache, _to_physical_domain(bs[i], cell_maps[i]), fh), 1:num_cells(trian))
  else
    cell_vals = lazy_map(i -> evaluate!(cache, bs[i], fh), 1:num_cells(trian))
  end
end
function _to_physical_domain(b::MomentBasedDofBasis, cell_map::Field)
  nodes = b.nodes
  f_moments = b.face_moments
  f_phys_nodes = cell_map(nodes)
  MomentBasedDofBasis(f_phys_nodes, f_moments, b.face_nodes)
end
function _to_physical_domain(b::LagrangianDofBasis, cell_map::Field)
  nodes = b.nodes
  f_phys_nodes = cell_map(nodes)
  LagrangianDofBasis(f_phys_nodes, b.dof_to_node, b.dof_to_comp, b.node_and_comp_to_dof)
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
