using Gridap
using Gridap.ReferenceFEs
using Gridap.Arrays
import Gridap.Arrays.evaluate!

function CF(V::FESpace, fh::FEFunction)
  trian = get_triangulation(V)
  b = get_fe_dof_basis(V)
  bs = get_data(b)
  cache = return_cache(testitem(bs), fh)
  cell_coords = get_cell_points(b).cell_phys_point
  cell_vals = lazy_map(i -> evaluate!(cache, bs[i], fh, cell_coords[i]),
                       1:num_cells(trian))
  CellField(V, cell_vals)
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
  vals = evaluate!(cf,field,b.nodes)
  ndofs = length(b.dof_to_node)
  T = eltype(vals)
  ncomps = num_components(T)
  ReferenceFEs._evaluate_lagr_dof!(c,vals,b.node_and_comp_to_dof,ndofs,ncomps)
end
