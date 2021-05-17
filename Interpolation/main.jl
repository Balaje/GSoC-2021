using Gridap
using GridapGmsh
using Gridap.ReferenceFEs
import Gridap: gradient

## This is from the 3D problem of iceberg stuff.
domain = (0, 1, 0, 1)
partition = (200,200)
model = CartesianDiscreteModel(domain,partition)
#model = GmshDiscreteModel("../testmesh.mesh")

T=Float64
order=1
pol=QUAD
reffe=LagrangianRefFE(T,pol,order)



# Functions u and f.
u(x) = sin(π*x[1])*sin(π*x[2])
f(x) = 2*π^2*sin(π*x[1])*sin(π*x[2])

## FE Spaces
#order = 1;
#reffe = ReferenceFE(lagrangian,Float64,1);
V0 = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary");
U = TrialFESpace(V0,0.);

# Degree and domain
degree = 2;
Ω = Triangulation(model);
dΩ = Measure(Ω,degree);

# Weak form
a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ;
b(v) = ∫( v*f )*dΩ;

op = AffineFEOperator(a,b,U,V0);

uh=solve(op)
