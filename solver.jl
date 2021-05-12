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
ω=2*π/60; # 40s incident wave.
N=5; # Modal expansion in the ocean
nev=20; #Number of eigenvalues
L=10000; #Shelf length
h=200; #Shelf thickness
d=0.9*h; #Submergence.
H=800; #Ocean depth
E=2.0e9;
ρᵢ=922.5;
ρₗ=1025;
ν=0.33;
EI=E*h^3/(12*(1-ν^2));
Lc=(EI/(ρₗ*9.8))^0.25;
tc=sqrt(ρₗ*Lc^6/(EI*H));
LL=L/Lc;
HH=H/Lc;
hh=h/Lc;
dd=d/Lc;
α=HH*(ω*tc)^2;
Ap=(9.8/(1im*ω));

# Solve the dispersion equation
k=dispersionEquations.dispersionfreesurface(α, N, HH);
k[1]=-k[1];
kd=dispersionEquations.dispersionfreesurface(α, N, HH-dd);
kd[1]=-kd[1];
print("Solved dispersion equations\n")

partition=(200,10);

# Build model for the ice-shelf.
#iceDomain=(0,L,-d,h-d);
#iceModel=CartesianDiscreteModel(iceDomain,partition);
#iceLabels=get_face_labeling(iceModel);
#Ωs=Triangulation(iceModel); #Build the triangulation
#add_tag_from_tags!(iceLabels,"neumannIce",[1,2,5])
#add_tag_from_tags!(iceLabels,"dirichletIce",[4,8])
#Γs₁=BoundaryTriangulation(iceModel,iceLabels,tags="neumannIce");
#nΓs₁=get_normal_vector(Γs₁);# Get the normal
#dΓs₁=Measure(Γs₁,2);
#reffe=ReferenceFE(lagrangian,VectorValue{2,Float64},1);
#Vsₕ=TestFESpace(iceModel,reffe,conformity=:H1,dirichlet_tags="dirichletIce",
#                dirichlet_masks=(true,true)); #Test space for Ice

## Build model for the cavity region
cavDomain=(0,LL,-HH,-dd);
cavModel=CartesianDiscreteModel(cavDomain,partition);
cavLabels=get_face_labeling(cavModel);
Ωf=Triangulation(cavModel); #Build the triangulation
add_tag_from_tags!(cavLabels,"neumannIce",[6]);
add_tag_from_tags!(cavLabels,"NonLocal",[3,7]);
Γf₃=BoundaryTriangulation(cavModel,cavLabels,tags="neumannIce"); # Shelf/Cavity
Γf₄=BoundaryTriangulation(cavModel,cavLabels,tags="NonLocal"); # Non-Local
reffe=ReferenceFE(lagrangian,Float64,1);
Vfₕ=TestFESpace(cavModel,reffe,conformity=:H1); #Test space for cavity.

# ------ Get the non-local boundary condition
Qϕ,χ=nonLocalBoundary.getMQχ(k, kd, HH, dd, N, Ap, cavModel, Γf₄, Vfₕ, Vfₕ);
print("Done computing non-local boundary condition\n")
# -----------------------------------------

# First attempt at solving the Elasticity Eigenvalue problem
# Question on interpolation?
#ξ,Vec=FEMSolvers.solveEigen(iceModel, Vsₕ, Vsₕ, 1, 10); # Returns the raw vectors.
#uh=FEFunction(Vsₕ, real(Vec[:,1]))
#writevtk(Ωs,"results",cellfields=["uh"=>uh]); #To visualize the solution.

# ------- Solve for the velocity potentials
#Diffraction Potential (with a rigid shelf)
K,f,op=FEMSolvers.getLaplaceMatEB(Ωf, Γf₃, Vfₕ, Vfₕ, Qϕ, χ, 0, LL, 0)
ϕ₀=K\f;
# Radiation potential from the Eigenmodes
μ=dispersionEquations.solveEigenEB(nev, LL);# Obtain the Eigenvalues for the beam equation
ϕₖ=zeros(ComplexF64,length(χ),nev)
for m=1:nev
    K,f,op=FEMSolvers.getLaplaceMatEB(Ωf, Γf₃, Vfₕ, Vfₕ, Qϕ, 0*χ, μ[m], LL, ω*Lc)
    ϕₖ[:,m]=K\f; # Using the linear algebra package (raw vector)
end
print("Done computing potentials\n")
# -----------------------------------------

# ------ Build and solve the reduced system
λ,K,B,AB,F=FEMSolvers.buildReducedSystem(μ, ϕ₀, ϕₖ, α, 1, dd, Γf₃, LL, ω, Vfₕ);
print("Done computing coefficients\n")
# ----------------------------------

## Construct the displacement
function u(x, μ, L, λ)
    nev=length(μ);
    ξ=0;
    for m=1:nev
        μₘ=μ[m];
        ηₘ = ((cos(L*μₘ) + cosh(L*μₘ))*(sin(μₘ*x) + sinh(μₘ*x))-
                (sin(L*μₘ) + sinh(L*μₘ))*(cos(μₘ*x) + cosh(μₘ*x)))/
                (cos(L*μₘ) + cosh(L*μₘ));
        ξ = ξ+λ[m]*ηₘ;
    end
    return ξ;
end

X=collect(range(0, LL, length=100));
U=zeros(ComplexF64,length(X),1);
for m=1:length(U)
    U[m]=u(X[m], μ, LL, λ);
end
# Construct the velocity potential
POT=ϕ₀+ϕₖ*λ;
uh=FEFunction(Vfₕ,real(POT[:,1]))
writevtk(Ωf,"results",cellfields=["uh"=>uh]); #To visualize the solution.
