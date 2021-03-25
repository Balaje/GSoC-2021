using Gridap
using GridapGmsh
import Gridap: gradient

## This is from the 3D problem of iceberg stuff.
domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)


# Functions u and f.
u(x) = x[1] + x[2]
f(x) = 0

## FE Spaces
order = 1;
reffe = ReferenceFE(lagrangian,Float64,order);
V0 = TestFESpace(model, reffe, conformity=:H1,dirichlet_tags="boundary");
U = TrialFESpace(V0,u);

# Degree and domain
degree = 2;
Ω = Triangulation(model);
dΩ = Measure(Ω,degree);

# Weak form
a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ;
b(v) = ∫( v*f )*dΩ;

op = AffineFEOperator(a,b,U,V0);
