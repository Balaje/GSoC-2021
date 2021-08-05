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
 interpolate_everywhere(Interpolatable(fₕ), V₂)  Works without these routines.
   Here fₕ ∈ V₁ ≠ V₂
 Check interpolable_1.jl

Functions for Optimizing interpolations:
- _cell_vals
- _to_physical_domain
"""

using Gridap.Helpers
using Gridap.Arrays

struct PushDofMap <: Map end
function Arrays.evaluate!(cache, ::PushDofMap, f::Dof, m::Field)
  @abstractmethod
end
function Arrays.evaluate!(cache, ::PushDofMap, f::AbstractArray{<:Dof}, m::Field)
  @abstractmethod
end

# Local implementations
function replace_nodes(f::LagrangianDofBasis, x)
  LagrangianDofBasis(x, f.dof_to_node, f.dof_to_comp, f.node_and_comp_to_dof)
end
function Arrays.return_cache(::PushDofMap, f::LagrangianDofBasis, m::Field)
  q = f.nodes
  return_cache(m, q)
end
function Arrays.evaluate!(cache, ::PushDofMap, f::LagrangianDofBasis, m::Field)
  q = f.nodes
  x = evaluate!(cache,m,q)
  replace_nodes(f,x)
end
function Arrays.lazy_map(::PushDofMap, cell_f::AbstractArray{<:LagrangianDofBasis}, cell_m::AbstractArray{<:Field})
  cell_q = lazy_map(f->f.nodes, cell_f)
  cell_x = lazy_map(evaluate, cell_m, cell_q)
  lazy_map(replace_nodes, cell_f, cell_x)
end

# CellData
function CellData.change_domain(a::CellDof, ::ReferenceDomain, ::PhysicalDomain)
  trian = get_triangulation(a)
  cell_m = get_cell_map(trian)
  cell_f_ref = get_data(a)
  cell_f_phys = lazy_map(PushDofMap(), cell_f_ref, cell_m)
  CellDof(cell_f_phys, trian, DomainStyle(a))
end
function FESpaces._cell_vals(V::SingleFieldFESpace, f::Interpolatable)
  fe_basis = get_fe_dof_basis(V)
  trian = get_triangulation(V)
  fh = f.uh
  s = CellData.change_domain(fe_basis, ReferenceDomain(), PhysicalDomain())
  bs = get_data(s)
  cache = return_cache(testitem(bs), fh)
  cell_vals = lazy_map(i-> evaluate!(cache, bs[i], fh), 1:num_cells(trian))
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
