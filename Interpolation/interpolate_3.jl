using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap
using Test

function Gridap.Arrays.return_cache(f::Gridap.Polynomials.QCurlGradMonomialBasis{D,T},x::Point) where {D,T}
  ndof = size(f.qgrad)[1]
  r = Gridap.Arrays.CachedArray(zeros(VectorValue{D,T},(ndof,)))
  xs = [x]
  cf = Gridap.Arrays.return_cache(f,xs)
  r, cf, xs
end
function Gridap.Arrays.evaluate!(cache,f::Gridap.Polynomials.QCurlGradMonomialBasis{D,T},x::Point) where {D,T}
  r, cf, xs = cache
  xs[1] = x
  v = Gridap.Arrays.evaluate!(cf,f,xs)
  ndof = size(v,2)
  Gridap.Arrays.setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end

p = QUAD
D = num_dims(QUAD)
et = Float64
f(x) = VectorValue([x[1]+x[2],0])
reffe = RaviartThomasRefFE(et, p, 0)

# Source RT Space.
model = CartesianDiscreteModel((0,1,0,1),(10,10))
V1 = FESpace(model, reffe)
fh = interpolate_everywhere(f, V1)

# Target RT Space
model = CartesianDiscreteModel((0,1,0,1),(10,10))
V2 = FESpace(model, reffe)

# Interpolate to the new space.
function get_rt_cell_dofs(fₕ, V₂)
  b = get_fe_dof_basis(V₂)
  fₕ_phys_point = x -> evaluate(fₕ, x)
  cp = get_cell_points(b)
  fₕ_phys_point_x = lazy_map(fₕ_phys_point, cp.cell_phys_point)
  lazy_map(x -> ReferenceFEs._eval_moment_dof_basis!(x, x, b.
    cell_dof[1]), fₕ_phys_point_x)
end

#gₕ = CellField(V2, get_rt_cell_dofs(fh, V2))


#@test evaluate(fh, Point(0.1,0.9)) ≈ evaluate(f, Point(0.1,0.9))
