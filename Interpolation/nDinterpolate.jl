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


xmin = 0
xmax = 1
ncells = 5


for D=1:3
    # Define domain and an affine function with arbitrary coeffs.
    domain = repeat([xmin,xmax], D)
    coeff0 = rand(Float64)
    coeffs = rand(SVector{D,Float64})
    f(x) = coeffs ⋅ SVector(Tuple(x)) + coeff0

    # Build the base model,
    base_partition = repeat([ncells], D)
    base_model = CartesianDiscreteModel(domain,base_partition)
    base_order = 1
    base_reffe = ReferenceFE(lagrangian, Float64, base_order)
    # Construct a FEfunction in the base_model.
    base_V1 = FESpace(base_model, base_reffe)
    global fh = interpolate_everywhere(f, base_V1) # Old FEFunbction
    writevtk(get_triangulation(fh),"source"*string(D),cellfields=["fh"=>fh])

    map_types = [:random, :sinusoidal]

    for map_type ∈ map_types
        # Second FE Space
        function rndm(p::Point, D=pt_to_rand, mt=map_type)
            vec = map(x -> collect(Tuple(x)), p)
            if(mt==:random)
                Random.seed!(get(D, p, 1234))
                per = -0.00625 + 0.0125*rand()
                new_vec = map(r -> (r≈1.0||r≈0.0) ? r : r+per, vec)
            elseif(mt==:sinusoidal)
                # s2π = Πᵢ[sin(2πxᵢ)]
                s2π = 0.1*prod(map(x -> sin(2π*x), vec)) # Appropriate scaling
                new_vec = map(x -> s2π+x, vec)
            end
            VectorValue(new_vec)
        end

        partition = map(x -> 4x, base_partition)
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
        order = 1
        reffe = ReferenceFE(lagrangian, Float64, order)
        V2 = FESpace(model, reffe)

        phys_point = get_cell_points(get_fe_dof_basis(V2)).cell_phys_point
        fh_phys_coords(x) = evaluate(fh, x)
        phys_point_fx = lazy_map(fh_phys_coords, phys_point)
        global gh = CellField( V2, phys_point_fx )

        writevtk(get_triangulation(gh),"target"*string(map_type)*string(D),cellfields=["fh"=>gh])

        #xs = [VectorValue(rand(D)) for i in 1:3]
        #for x in xs
        #@show D
        #fx = f(x)
        #fhx = evaluate(fh, x)
        #ghx = evaluate(gh, x)
        #end
    end
end
