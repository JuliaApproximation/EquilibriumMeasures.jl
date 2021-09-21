using EquilibriumMeasures, StaticArrays, ClassicalOrthogonalPolynomials, FillArrays, BlockArrays, LazyBandedMatrices, LinearAlgebra, IntervalSets, Test
import EquilibriumMeasures: EquilibriumMeasureMoment
using ForwardDiff: jacobian
using DomainSets: components

@testset "one-interval" begin
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
    μ, b = equilibriummeasure(x -> (x-3)*(x-2)*(1+x)*(2+x)*(3+x)*(2x-1)/20; a=SVector(-3,-2), returnendpoint=true)
    @test sum(μ) ≈ 1 atol=1E-13
    
    # Deflate found solution and find another one whilst damping the Newton step
    μ = equilibriummeasure(x -> (x-3)*(x-2)*(1+x)*(2+x)*(3+x)*(2x-1)/20; a=SVector(-3,-2), knownsolutions=[b], dampening=0.3)
    @test sum(μ) ≈ 1 atol=1E-6
end

@testset "two-interval" begin
    # see EquilibriumMeasureExamples.nb
    V = x -> 0.7*(x^4 - 2x^3 - x^2 + 2x)
    a_ex = SVector(-1.0637226766068189, 0.2659671162729329, 0.7340328837270674, 2.063722676606819)
    @test norm(EquilibriumMeasureMoment(V)(a_ex)) ≤ 1E-13

    @time μ = equilibriummeasure(V; a=SVector(-1,0.25,0.7,2))
    @test all(components(axes(μ,1).domain) .≈ (a_ex[1]..a_ex[2], a_ex[3]..a_ex[4]))

    
    V = x -> (x-3)*(x-2)*(1+x)*(2+x)*(3+x)*(2x-1)/150
    a_ex = SVector(-3.0782937003244046, 0.9082915776561011, 1.8384672185745663, 3.0223539638196737)
    @test norm(EquilibriumMeasureMoment(V)(a_ex)) ≤ 1E-13
    
    @time μ = equilibriummeasure(V; a=SVector(-3,1,2,3))
    @test vcat(leftendpoint.(components(axes(μ,1).domain))...) ≈ a_ex[[1,3]]
    @test vcat(rightendpoint.(components(axes(μ,1).domain))...) ≈ a_ex[[2,4]]
end