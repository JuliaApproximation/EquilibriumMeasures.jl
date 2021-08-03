using EquilibriumMeasures, Test

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
end
    


