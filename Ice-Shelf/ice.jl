using Gridap
import Gridap: ∇
using Arpack
using Roots
using SparseArrays

# Add packages
include("dispersion.jl");

include("nonLocal.jl");

include("FEMSolve.jl");
