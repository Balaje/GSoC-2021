# Gridap GSoC 2021

Contains some code related to my Google Summer of Code. I mostly fill
this up with experiments and to test concepts. Checkout `DETAILS.md`
for more details and also the blog post at the end.

## Evaluating function on an arbitrary point.

Extend `evaluate` to handle arbitrary points. Build `kdtree` using the coordinates and search the tree for the cell containing the point. Already under development in Gridap. It works:

```julia
reffe = ReferenceFE(lagrangian, Float64, 1)
model = CartesianDiscreteModel((0,1,0,1), (20,20))
V1 = FESpace(model, reffe)
fh = interpolate_everywhere(f, V1)

evaluate(fh, Point(rand(2))
fh(Point(rand(2)) # Also works!
```

## Use `evaluate` to interpolate between two meshes.

Interpolation can be done by using `evaluate` from before. The simplest one is 

```julia
phys_point = get_cell_points(get_fe_dof_basis(V2)).cell_phys_point
fh_phys_coords(x) = evaluate(fh, x)
phys_point_fx = lazy_map(fh_phys_coords, phys_point)
```

| ![2d](Images/2d.png) | ![3d](Images/3d.png) |
| -- | -- |

One of my mentors came up with this neat solution.
```julia
model = CartesianDiscreteModel((0,1,0,1), (40,40))
V2 = FESpace(model, reffe)
gh = interpolate_everywhere(x->fh(x), V2)
```
But this is a little bit slow since it builds the trees everytime a new point is encountered. There are a couple of ways to do it. Most define a method `_cell_vals(b::CellDof, fh)` where `b` is the cell basis of the target space `V2` obtained using `b = get_fe_dof_basis(V2)`. Then evaluating them in the new cell  points

```julia
function _cell_vals(b::CellDof, fh::CellField)
  # Compute cell dof values of fh on the new mesh
  return cell_vals # cell-wise dof values
end
```

This can be done in many ways and `README.md` will be soon updated once the optimal way is implemented. 

## Blog posts

[Main Page](https://balaje.github.io/gsoc/index.html)

1. [Google Summer of Code with Gridap.jl](https://balaje.github.io/2021/05/19/Google-Summer-Code.html)
2. [Google Summer of Code 2021 - Community Bonding](https://balaje.github.io/2021/06/05/GSoC-Week-0.html)
3. [Interpolation for Lagrangian Elements - I](https://balaje.github.io/2021/07/05/GSoC-Week-1.html)
4. [Interpolation for Lagrangian Elements - II](https://balaje.github.io/2021/07/05/GSoC-Week-2.html)
