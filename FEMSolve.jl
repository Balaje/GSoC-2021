module FEMSolvers
using Gridap
using Arpack

# Solve the eigenvalue problem. [ε is a function to store the symmetric gradient]
function solveEigen(iceModel,Vh,Vh0,order,N)
    λ=3500;
    ν=0.3;
    μ=λ*(1-2*ν)/(2*ν);
    ρᵢ=1;
    # The stress tensor
    σ(ε)=λ*tr(ε)*one(ε) + 2*μ*ε;
    Ω=Triangulation(iceModel);
    dΩ=Measure(Ω,2*order);
    # Weak form
    a(u,v) = ∫(ε(v)⊙(σ∘ε(u)))*dΩ;
    m(u,v) = ∫(ρᵢ*u⋅v)*dΩ;
    l(v) = 0;
    # Get the matrices
    opK = AffineFEOperator(a,l,Vh,Vh0);
    opM = AffineFEOperator(m,l,Vh,Vh0);

    K=opK.op.matrix;
    M=opM.op.matrix;

    ξ,Vec = eigs(K,M; nev=N, which=:SM);
    return ξ,Vec;
end

end
