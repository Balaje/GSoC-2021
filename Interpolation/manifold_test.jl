using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Arrays
using Gridap
using Test
using StaticArrays

@testset "evaluating functions" begin
  for D in 1:3
    xmin = 0
    xmax = 1
    domain = repeat([xmin, xmax], D)
    ncells = 3
    partition = repeat([ncells], D)
    base_model = CartesianDiscreteModel(domain, partition)
    order = 2
    reffe = ReferenceFE(lagrangian, Float64, order)

    for M ∈ [:hypercubes, :simplices]
      model = (M == :hypercubes) ? base_model : simplexify(base_model)

      V = FESpace(model, reffe)

      coeff0 = rand(Float64)
      coeffs = rand(SVector{D,Float64})
      f(x) = coeffs ⋅ SVector(Tuple(x)) + coeff0
      # TODO: use this mechanism instead to project
      # Francesc Verdugo @fverdugo 13:11
      # a(u,v) = ∫( u*v )dΩ
      # l(v) = a(f,v)
      # Solve a fe problem with this weak form
      # See also tutorial 10, "Isotropic damage model", section "L2
      # projection", function "project"
      fh = interpolate_everywhere(f, V)
      fhcache = return_cache(fh, VectorValue(zeros(D)...))

      # Test Random points
      xs = [VectorValue(rand(D)...) for i in 1:10]
      for x in xs
        x = VectorValue(rand(D)...)
        fx = f(x)
        fhx = evaluate!(fhcache, fh, x)
        @test fhx ≈ fx
      end
      fhxs = fh(xs)
      @test fhxs ≈ f.(xs)

      nv = num_vertices(model) # Number of vertices
      nf = num_faces(model,D-1) # Number of faces
      trian = Triangulation(model)
      topo = GridTopology(model)

      pts = get_vertex_coordinates(topo) # Vertex coordinates
      face_nodes = get_face_nodes(model, D-1) # face-to-node numbering
      face_coords = lazy_map(Broadcasting(Reindex(pts)), face_nodes) # Get LazyArray of coordinates of face

      # Test a random vertex from the triangulation
      pt = pts[rand(1:nv)]
      fhpt = evaluate!(fhcache, fh, pt)
      @test fhpt .≈ f(pt)

      # Test a random point lying on the face of the polytope
      face_coord = face_coords[rand(1:nf)]
      λ = rand(length(face_coord));
      λ = (D > 1) ? λ./sum(λ) : λ
      pt = face_coord ⋅ λ # Point on the face
      fhpt = evaluate!(fhcache, fh, pt)
      @test fhpt .≈ f(pt)

      # Test with CellPoint
      # Build cell_to_fxs manually
      cache1,cache2 = fhcache
      ncells = num_cells(model)
      x_to_cell(x) = CellData._point_to_cell!(cache1, x)
      point_to_cell = map(x_to_cell, xs)
      cell_to_points, point_to_lpoint = make_inverse_table(point_to_cell, ncells)
      cell_to_xs = lazy_map(Broadcasting(Reindex(xs)), cell_to_points)
      cell_to_f = get_array(fh)
      cell_to_fxs = lazy_map(evaluate, cell_to_f, cell_to_xs)

      # Now build CellPoint with xs instead of building cell_to_xs
      cell_point_xs = compute_cell_points_from_vector_of_points(xs, trian, PhysicalDomain())
      cell_point_fxs = evaluate(fh, cell_point_xs)
      @test cell_point_fxs ≈ cell_to_fxs

    end
  end
end

@testset "evaluating functions on manifolds" begin
  # TODO: Test for D=2
  for D in 2:4
    # Create domain
    xmin=0
    xmax=1

    domain = repeat([xmin,xmax], D)
    ncells = 3
    partition = repeat([ncells], D)

    order = 2
    reffe = ReferenceFE(lagrangian, Float64, order)

    model = CartesianDiscreteModel(domain, partition)
    bgmodel = simplexify(model)
    labels = get_face_labeling(bgmodel)
    bgface_to_mask = get_face_mask(labels, "boundary", D-1) # Get the boundary of the bgmodel
    bmodel = BoundaryDiscreteModel(Polytope{D-1}, bgmodel, bgface_to_mask)

    coeff0 = rand(Float64)
    coeffs = rand(SVector{D, Float64})
    f(x) = coeffs ⋅ SVector(Tuple(x)) + coeff0
    V = FESpace(bmodel, reffe)
    fh = interpolate_everywhere(f, V)

    xs = [VectorValue(rand(D-1)..., xmax) for i in 1:10]
    ys = [VectorValue(xmin, rand(D-1)...) for i in 1:10]
    for i in 1:length(xs)
      x = xs[i]
      fx = f(x)
      fhx = evaluate(fh, x)
      @test fx ≈ fhx

      y = ys[i]
      fy = f(y)
      fhy = evaluate(fh, y)
      @test fy ≈ fhy
    end

    fhxs = fh.(xs)
    fhys = fh.(ys)
    @test fhxs ≈ f.(xs)
    @test fhys ≈ f.(ys)

    # Random points outside bmodel, but inside model
    pt = [VectorValue(rand(D)) for i in 1:5]
    try
      fh.(pt)
    catch
      @test 1==1
    end

  end
end
