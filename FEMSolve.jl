module FEMSolvers
using Gridap
import Gridap: ∇
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
    # Use Arpack to solve the Eigenvalue problem
    ξ,Vec = eigs(K,M; nev=N, which=:SM);
    return ξ,Vec;
end


# Solve the velocity potentials.
function getLaplaceMatEB(Ω, Γ₃, Vh, Vh0, QΦ, χ, μₘ, L, ω)
    # Measure of the domains
    dΩ=Measure(Ω,2);
    dΓ₃=Measure(Γ₃,2); #Interface boundary.

    η(x) = ω*((cos(L*μₘ) + cosh(L*μₘ))*(sin(μₘ*x[1]) + sinh(μₘ*x[1]))-
            (sin(L*μₘ) + sinh(L*μₘ))*(cos(μₘ*x[1]) + cosh(μₘ*x[1])))/
            (cos(L*μₘ) + cosh(L*μₘ));

    a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ;
    b(v) = ∫(η*v)*dΓ₃;
    op=AffineFEOperator(a,b,Vh,Vh0);
    K=op.op.matrix+QΦ;
    f=op.op.vector+(ω==0)*χ[1,:];
    return K,f,op
end

end
