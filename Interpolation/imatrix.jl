#include("interpolable.jl")

include("interpolable_include.jl")
using SparseArrays

p = TRI
D = num_dims(QUAD)
et = Float64
source_model = simplexify(CartesianDiscreteModel((0,1,0,1),(10,10)))
f(x) = x[1] + x[2]
reffe = LagrangianRefFE(et, p, 2)
V₁ = FESpace(source_model, reffe, conformity=:H1)
fh = interpolate_everywhere(f, V₁)
# Target Lagrangian Space
reffe = LagrangianRefFE(et, p, 2)
model = simplexify(CartesianDiscreteModel((0,1,0,1),(20,20)))
V₂ = FESpace(model, reffe, conformity=:H1)


function interpmatrix(V₂, V₁)
  field = get_fe_basis(V₁)
  fields = get_array(field)
  σ₁ = get_cell_dof_ids(V₁)
  M = num_free_dofs(V₁)

  b = get_fe_dof_basis(V₂)
  b = change_domain(b, ReferenceDomain(), PhysicalDomain())
  bs = get_data(b)
  σ₂ = get_cell_dof_ids(V₂)
  N = num_free_dofs(V₂)

  cache = CellData._point_to_cell_cache(get_triangulation(V₁));

  A = spzeros(M,N)
  cell_nodes = lazy_map(x->x.nodes,bs)
  cell_basis_vals = lazy_map(node -> lazy_map(pt->evaluate(field,pt), node), cell_nodes)
  cell_inds = lazy_map(node -> lazy_map(pt->CellData._point_to_cell!(cache,pt), node), cell_nodes);

  for m ∈ 1:num_cells(get_triangulation(V₂))
    for i ∈ 1:length(σ₁[cell_inds[m]])
      A[σ₁[cell_inds[m]][i],σ₂[m][i]] = cell_basis_vals[m][i]
    end
  end
  A
end
