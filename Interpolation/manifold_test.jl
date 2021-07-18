using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Arrays
using Gridap
using Test
using StaticArrays


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
    pt = [VectorValue(rand(D)) for i in 1:3]
    try
      fh.(pt)
    catch
      @test 1==1
    end

  end
end
