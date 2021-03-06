using Gridap
using Gridap.Algebra
using Gridap.CellData
using Gridap.Fields
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.MultiField
using Gridap.ReferenceFEs
using Gridap.Visualization
using StaticArrays
using Random
using Test


xmin = 0
xmax = 1
ncells = 5


@testset begin
  for D=1:3
    # Define domain and an affine function with arbitrary coeffs.
    domain = repeat([xmin,xmax], D)
    coeff0 = rand(Float64)
    coeffs = rand(SVector{D,Float64})
    f(x) = coeffs ⋅ SVector(Tuple(x)) + coeff0

    # Build the base model,
    base_partition = repeat([ncells], D)
    base_model = CartesianDiscreteModel(domain,base_partition)
    base_order = 2
    base_reffe = ReferenceFE(lagrangian, Float64, base_order)
    # Construct a FEfunction in the base_model.
    base_V1 = FESpace(base_model, base_reffe)
    fh = interpolate_everywhere(f, base_V1) # Old FEFunbction
    #writevtk(get_triangulation(fh),"source"*string(D),cellfields=["fh"=>fh])

    map_types = [:random, :sinusoidal]

    for map_type ∈ map_types
      # Second FE Space
      function rndm(p::Point, D=pt_to_rand, mt=map_type)
        vec = map(x -> collect(Tuple(x)), p)
        if(mt==:random)
          Random.seed!(get(D, p, 1234))
          per = -1. /(8. *ncells) + 2. /(8. *ncells)*rand()
          new_vec = map(r -> (r≈1.0||r≈0.0) ? r : r+per, vec)
        elseif(mt==:sinusoidal)
          # s2π = Πᵢ[sin(2πxᵢ)]
          s2π = 0.02*prod(map(x -> sin(2π*x), vec)) # Appropriate scaling
          new_vec = map(x -> s2π+x, vec)
        end
        VectorValue(new_vec)
      end

      partition = base_partition
      if(map_type==:random)
        # Construct a uniform partition first to build the perturbation
        model_uni = CartesianDiscreteModel(domain,partition)
        # Random permutation of indices
        pt_to_ind = randperm(num_nodes(model_uni))
        # Build a dictionary with random indices for points
        pt_to_rand = Dict(zip(get_node_coordinates(model_uni),pt_to_ind))
        model = CartesianDiscreteModel(domain, partition; map=rndm)
      else
        pt_to_rand = Nothing # Create empty dict for sinusoidal maps
        model = CartesianDiscreteModel(domain, partition; map=rndm)
      end

      # Build the new FESpace
      order = 2
      reffe = ReferenceFE(lagrangian, Float64, order)
      V2 = FESpace(model, reffe)

      gh = interpolate_everywhere_non_compatible_trian(fh, V2)

      #writevtk(get_triangulation(gh),"target"*string(map_type)*string(D),cellfields=["fh"=>gh])

      xs = [VectorValue(rand(D)) for i in 1:3]
      for x in xs
        fx = f(x)
        fhx = evaluate(fh, x)
        ghx = evaluate(gh, x)
        @test fx ≈ fhx ≈ ghx
      end
    end
  end
end
