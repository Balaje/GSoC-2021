using Gridap
using GridapGmsh

include("dispersion.jl");
include("nonLocal.jl");

using .nonLocalBoundary
using .dispersionEquations

# Some parameters
H=5; #Ocean depth
omega=2*π/100; # 100s incident wave.
alph=omega^2;
N=5; # Modal expansion in the ocean
Ap=(10/1im*omega);
L=10; #Shelf length
h=1; #Shelf thickness
d=0.9*h; #Submergence.

# Solve the dispersion equation
k=dispersionEquations.dispersionfreesurface(alph, N, H);
kd=dispersionEquations.dispersionfreesurface(alph, N, H-d);

# Build model for the cavity and ice-shelf.
cavDomain=(0,L,-H,-d);
iceDomain=(0,L,-d,h-d);
partition=(100,100);
iceModel=CartesianDiscreteModel(iceDomain,partition);
cavModel=CartesianDiscreteModel(cavDomain,partition);

## Label the boundaries for the cavity region
iceLabels=get_face_labeling(iceModel);
cavLabels=get_face_labeling(cavModel);

# --Set labels and triangulation for the cavity--
Ωf=Triangulation(cavModel);
add_tag_from_tags!(cavLabels,"neumannIce",[6]);
add_tag_from_tags!(cavLabels,"NonLocal",[3,7]);
Γfₙ=BoundaryTriangulation(cavModel,cavLabels,tags="neumannIce");
Γfₙₗ=BoundaryTriangulation(cavModel,cavLabels,tags="NonLocal");

# --Set labels and triangulation for the ice--
Ωs=Triangulation(iceModel); #Build the triangulation
add_tag_from_tags!(iceLabels,"neumannIce",[1,2,5])
add_tag_from_tags!(iceLabels,"dirichletIce",[4,8])
Γsₙ=BoundaryTriangulation(iceModel,iceLabels,tags="neumannIce");

#Build the spaces for the ice and cavity
reffe=ReferenceFE(lagrangian,VectorValue{2,Float64},1);
Vsₕ=TestFESpace(iceModel,reffe,conformity=:H1,dirichlet_tags="dirichletIce",dirichlet_masks=(true,true)); #Test space for Ice
reffe=ReferenceFE(lagrangian,Float64,1);
Vfₕ=TestFESpace(cavModel,reffe,conformity=:H1); #Test space for cavity.






#A,M,f,g=nonLocalBoundary.getMAT(k, kd, H, d, N, Ap);
