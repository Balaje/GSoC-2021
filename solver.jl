using GridapGmsh
using Gridap: ∇
using Gridap
using SparseArrays
using Roots

include("dispersion.jl");
include("nonLocal.jl");
include("FEMSolve.jl");

using .nonLocalBoundary
using .dispersionEquations
using .FEMSolvers

# Some parameters
H=5; #Ocean depth
ω=2*π/40; # 100s incident wave.
α=ω^2;
N=5; # Modal expansion in the ocean
Ap=(10/1im*ω);
L=20; #Shelf length
h=1; #Shelf thickness
d=0.9*h; #Submergence.

# Solve the dispersion equation
k=dispersionEquations.dispersionfreesurface(α, N, H);
kd=dispersionEquations.dispersionfreesurface(α, N, H-d);

partition=(100,10);

# Build model for the ice-shelf.
iceDomain=(0,L,-d,h-d);
iceModel=CartesianDiscreteModel(iceDomain,partition);
iceLabels=get_face_labeling(iceModel);
Ωs=Triangulation(iceModel); #Build the triangulation
add_tag_from_tags!(iceLabels,"neumannIce",[1,2,5])
add_tag_from_tags!(iceLabels,"dirichletIce",[4,8])
Γs₁=BoundaryTriangulation(iceModel,iceLabels,tags="neumannIce");
nΓs₁=get_normal_vector(Γs₁);# Get the normal
dΓs₁=Measure(Γs₁,2);
reffe=ReferenceFE(lagrangian,VectorValue{2,Float64},1);
Vsₕ=TestFESpace(iceModel,reffe,conformity=:H1,dirichlet_tags="dirichletIce",
                dirichlet_masks=(true,true)); #Test space for Ice

## Build model for the cavity region
cavDomain=(0,L,-H,-d);
cavModel=CartesianDiscreteModel(cavDomain,partition);
cavLabels=get_face_labeling(cavModel);
Ωf=Triangulation(cavModel); #Build the triangulation
add_tag_from_tags!(cavLabels,"neumannIce",[6]);
add_tag_from_tags!(cavLabels,"NonLocal",[3,7]);
# *) Boundary
Γf₃=BoundaryTriangulation(cavModel,cavLabels,tags="neumannIce");
nΓf₃=get_normal_vector(Γf₃);# Get the normal
Γf₄=BoundaryTriangulation(cavModel,cavLabels,tags="NonLocal");
reffe=ReferenceFE(lagrangian,Float64,1);
Vfₕ=TestFESpace(cavModel,reffe,conformity=:H1); #Test space for cavity.

# ------ Get the non-local boundary condition
Qϕ,χ=nonLocalBoundary.getMQχ(k, kd, H, d, N, Ap, cavModel,Γf₄,Vfₕ,Vfₕ);
# -----------------------------------------

# First attempt at solving the Elasticity Eigenvalue problem
# Question on interpolation?
#ξ,Vec=FEMSolvers.solveEigen(iceModel, Vsₕ, Vsₕ, 1, 10); # Returns the raw vectors.
#uh=FEFunction(Vsₕ, real(Vec[:,1]))
#writevtk(Ωs,"results",cellfields=["uh"=>uh]); #To visualize the solution.

# ------- Solve for the velocity potentials
#Diffraction Potential
K,f,op=FEMSolvers.getLaplaceMatEB(Ωf, Γf₃, Vfₕ, Vfₕ, Qϕ, χ, μ[1], L, 0)
f*=(-1im)
ϕ₀=K\f;
# Radiation potential from the Eigenmodes
nev=10; #Number of eigenvalues
μ=dispersionEquations.solveEigenEB(nev, L);# Obtain the Eigenvalues for the beam equation
ϕₖ=zeros(ComplexF64,length(χ),nev)
for m=1:nev
    K,f,op=FEMSolvers.getLaplaceMatEB(Ωf, Γf₃, Vfₕ, Vfₕ, Qϕ, χ, μ[m], L, ω)
    f*=(-1im)
    ϕₖ[:,m]=K\f; # Using the linear algebra package (raw vector)
    ϕₕ=FEFunction(Vfₕ,imag(ϕₖ[:,m]));
end
writevtk(Ωf,"results",cellfields=["uh"=>ϕₕ]); #To visualize the solution.
# -----------------------------------------
