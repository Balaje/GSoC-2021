# Program to generate the reflection coefficient plot in the real frequency space
L = 10000.
H = 800.
h = 200.
nev = 10
N = 5

ω = 2π/100;
Hmat, F, Ref, RefModes, RefDiff, X, U, Lc = solveIceVibration(10000,200,800,10,5,ω);
plotIce(X,U,ω,[-2.2,2.2],Ref)


# Get the coefficients on the new space
H_new, F_new, λ_new, ω_new = coarse2fine(2π/400, 2π/20, 10 ,200, L, h, H, nev, N)
plotMode(ω_new, λ_new, 4)
