using EquilibriumMeasures, StaticArrays, ClassicalOrthogonalPolynomials, Test

@testset "EquilibriumMeasures" begin
    μ = equilibriummeasure(x -> x^2)
    @test first(axes(μ,1)) ≈ -sqrt(2) atol=1E-13
    @test last(axes(μ,1)) ≈ sqrt(2) atol=1E-13
    @test sum(μ) ≈ 1 atol=1E-13
    @test μ[0.1] ≈ sqrt(2-0.1^2)/π atol=1E-13

    μ = equilibriummeasure(x -> x^4)
    @test first(axes(μ,1)) ≈ -1.0745699318235422 atol=1E-13
    @test last(axes(μ,1)) ≈ 1.0745699318235422 atol=1E-13
    @test sum(μ) ≈ 1 atol=1E-13
    @test μ[0.1] ≈ 0.40005825716887744

    μ = equilibriummeasure(x -> x^4 + x^3 - x^2)
    @test sum(μ) ≈ 1 atol=1E-13
    @test μ[0.1] ≈ 0.19674192408704289

    μ = equilibriummeasure(x -> (x-3)*(x-2)*(1+x)*(2+x)*(3+x)*(2x-1)/250; a=SVector(-3,3))
    @test sum(μ) ≈ 1 atol=1E-13

    # non-unique minimiser
    μ = equilibriummeasure(x -> (x-3)*(x-2)*(1+x)*(2+x)*(3+x)*(2x-1)/20; a=SVector(-3,-2))
    @test sum(μ) ≈ 1 atol=1E-13
end

@testset "two-interval" begin
    V = x -> x^4 - 10x^2
    xx = range(-4,4; length=1000)
    plot(xx, V.(xx))
    a,b,c,d = -3,-1,1,3
    W = PiecewiseInterlace(Weighted(chebyshevu(a..b)), Weighted(chebyshevu(c..d)))
    T = PiecewiseInterlace(chebyshevt(a..b), chebyshevt(c..d))
    x = axes(W,1)
    H = T \ inv.(x .- x') * W
    # H[Block.(2:
    E1 = PseudoBlockArray([Eye(2); Zeros(∞,2)], (axes(W,2), Base.OneTo(2)))
    H2 = BlockHcat(E1, H)
    H2 \ (T \ V.(x))
end