using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Arrays
using Gridap
using Test

p = QUAD
D = num_dims(QUAD)
et = Float64
f(x) = VectorValue([x[1], x[2]])
reffe = RaviartThomasRefFE(et, p, 0)

# Source RT Space.
model = CartesianDiscreteModel((0,1,0,1),(10,10))
V1 = FESpace(model, reffe, conformity=:HDiv)
fh = interpolate_everywhere(f, V1)

# Target RT Space
reffe = RaviartThomasRefFE(et, p, 1)
model = CartesianDiscreteModel((0,1,0,1),(10,10))
V2 = FESpace(model, reffe, conformity=:HDiv)

#### RT-INTERPOLATION
function rt_interpolate(f, V2)
  b = get_fe_dof_basis(V2)
  trian = get_triangulation(V2)
  cell_maps = get_cell_map(trian)
  phys_point = get_cell_points(b).cell_phys_point
  phys_point_fx = lazy_map(x -> evaluate(f, x), phys_point)
  bs = get_data(b)

  cell_dof_vals = lazy_map(i -> _eval_moment_dof_basis(phys_point_fx[i], bs[i]), 1:num_cells(trian))

  FEFunction(V2, gather_free_values(V2,cell_dof_vals))
end

function _eval_moment_dof_basis(vals::AbstractVector, b)
  T = eltype(vals)
  ndofs = num_dofs(b)
  dofs = zeros(eltype(T),ndofs)
  ReferenceFEs._eval_moment_dof_basis!(dofs, vals, b)
  dofs
end

# Now to test it
gh = rt_interpolate(fh, V2);
pt = Point(rand(2))
@test evaluate(fh, pt) ≈ evaluate(gh, pt)
