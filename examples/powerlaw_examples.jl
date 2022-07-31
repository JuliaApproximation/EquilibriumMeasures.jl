# This file acts as a simple temporary companion to the paper at https://arxiv.org/abs/2109.00843
# until the full functionality described in the paper is more properly integrated into EquilibriumMeasures.jl.

# Import required packages
using ClassicalOrthogonalPolynomials, QuasiArrays, ContinuumArrays, BandedMatrices, LazyArrays, LinearAlgebra,
        FastTransforms, ArrayLayouts, Test, FillArrays, BlockArrays, LazyBandedMatrices, ForwardDiff, HypergeometricFunctions, SpecialFunctions
import ClassicalOrthogonalPolynomials: Clenshaw, recurrencecoefficients, clenshaw, paddeddata, jacobimatrix, oneto, Weighted, MappedOPLayout
import LazyArrays: ApplyStyle
import QuasiArrays: MulQuasiMatrix
import Base: OneTo
import ContinuumArrays: MappedWeightedBasisLayout, Map, WeightedBasisLayout

# These definitions allow the use of radial Jacobi bases
struct QuadraticMap{T} <: Map{T} end
struct InvQuadraticMap{T} <: Map{T} end
QuadraticMap() = QuadraticMap{Float64}()
InvQuadraticMap() = InvQuadraticMap{Float64}()
Base.getindex(::QuadraticMap, r::Number) = 2r^2-1
Base.axes(::QuadraticMap{T}) where T = (Inclusion(0..1),)
Base.axes(::InvQuadraticMap{T}) where T = (Inclusion(-1..1),)
Base.getindex(d::InvQuadraticMap, x::Number) = sqrt((x+1)/2)
ContinuumArrays.invmap(::QuadraticMap{T}) where T = InvQuadraticMap{T}()
ContinuumArrays.invmap(::InvQuadraticMap{T}) where T = QuadraticMap{T}()
Base.getindex(d::QuadraticMap, x::Inclusion) = d
Base.getindex(d::InvQuadraticMap, x::Inclusion) = d

# Test the radial Jacobi basis
@testset "Quadratic map" begin
    P = Jacobi(0.4,0.1)[QuadraticMap(),:]
    Pn = Jacobi(0.4,0.1)[QuadraticMap(),1:10]
    Pn2 = P[:,Base.OneTo(10)]
    @test MemoryLayout(P) == MemoryLayout(Pn) == MemoryLayout(Pn2) == ContinuumArrays.MappedBasisLayout()
    x = axes(P,1)
    un = Pn \ (2 * x .^2 .- 1)
    un2 = Pn2 \ (2 * x .^2 .- 1)
    u = P \ (2 * x .^2 .- 1)
    @test un ≈ un2 ≈ u[1:10]
end

# Define the coefficients for the recurrence relationship
function frakca(α,β,d,m,n,z)
    return -(((2*m+4*n-α)*(d-2*(1+m+n)+α)*(-8*n^2*I-4*z+4*m^2*z+16*n^2*z+4*n*α*I-8*n*z*α+z*α^2-2*α*β*I+2*β^2*I+4*m*(-2*n*I+4*n*z-z*α+β*I)+d*I*(2+2*m-α+2*β)))/(2*(1+n)*(-2+2*m+4*n-α)*(2+2*m+2*n-α+β)*(d+2*m+2*n-α+β)))
end
function frakcc(α,β,d,m,n)
    return (((2+2*m+4*n-α)*(d-2*(m+n)+α)*(d-2*(1+m+n)+α)*(-2+2*n-β)*(d-2*n+β))/(4*n*(1+n)*(-2+2*m+4*n-α)*(2+2*m+2*n-α+β)*(d+2*m+2*n-α+β)))
end
# For simplicity we use the dense recurrence instead of the truncated-banded version but the truncated-banded version is to be preferred when computational efficiency is a concern.
function densehyper2(P, x, m, n, α, β, d)
    upper = π^(d/2)*gamma(β/2+1)*gamma((β+d)/2)*gamma(m+n-(α+d)/2+1)
    lower = gamma(d/2)*gamma(n+1)*gamma(β/2-n+1)*gamma((β-α)/2+m+n+1) 
    return (P \ ((upper/lower).*HypergeometricFunctions._₂F₁general2.(n-β/2,-m-n+(α-β)/2,d/2,abs.(x.^2))))
end
function ContiguousRecurrenceDense(α::Real, β::Real, d, a, b, m, NN::Integer)
    buffer = 30
    NN = NN+buffer
    # the space of RHS expansion
    P = Jacobi(a,b)[QuadraticMap(),1:NN+buffer];
    x = axes(P,1)
    # initializing the operator and dummy container
    Q = BandedMatrix{eltype(P)}(undef,(NN,NN),(NN,NN));
    # define jacobimatrix
    J = ((jacobimatrix(Jacobi(a,b))+I))[1:NN,1:NN]./2
    # generate the entries iteratively
    Q[1:NN,1] = densehyper2(P,x,l,0,α,β,d)[1:NN]
    Q[1:NN,2] = densehyper2(P,x,l,1,α,β,d)[1:NN]
    for n = 2:NN-1
        Q[1:NN,n+1] = frakca(α,β,d,m,n-1,J)*Q[1:NN,n]+frakcc(α,β,d,m,n-1)*Q[1:NN,n-1]
    end
    return Q[1:NN-buffer,1:NN-buffer]
end

# This defines an example problem for which we find the equilibrium measure
α = 1.31 # define attractive power
β = 1.23 # define repulsive power
l = 1 # set appropriate basis shift parameter
M = 1 # asks for probability measures
d = 2 # dimension of the problem, this is a 2D disk example

# This defines the basis of Jacobi polynomials.
# We choose the basis in which the attractive operator is banded
a = l-(β+d)/2 # Jacobi polynomial basis parameter
b = (d-2)/2   # Jacobi polynomial basis parameter
P = Jacobi(a,b)[QuadraticMap(),:]
x = axes(P,1)

# Compute the attractive and repulsive operators
op1 = ContiguousRecurrenceDense(α,α,d,a,b,l,100)
op2 = ContiguousRecurrenceDense(α,β,d,a,b,l,100)

# Energy computing helper function
function djacobiintegralnormalization(a,d,ρ)
    return π^(d/2)*gamma(a+1)/gamma(a+d/2+1)*ρ[1]
end
function computeEnerg(R, op1, op2, d, α, β, a, b, M)
    # Combine operators
    op = (R^(α)/α*op1-R^(β)/β*op2) # note: this is the non-regularized version
    # initial energy
    E = zeros(eltype(op1),size(op1,1))
    E[1] = one(eltype(op1))
    energ = qr(op) \ E
    # normalize
    ρ = (M/(djacobiintegralnormalization(a,d,energ))).*energ
    return (op*ρ)[1]
end
enplothelp(x) = computeEnerg(x,op1,op2,d,α,β,a,b,M)

# Measure computing helper function
function computemeas(R,op1,op2,d,α,β,a,b,M)
    # define operator
    op = R^(α)/α*op1-R^(β)/β*op2
    # initial energy
    E = zeros(eltype(op1),size(op1,1))
    E[1] = one(eltype(op1))
    energ = qr(op) \ E
    # normalize
    ρ = (M/(djacobiintegralnormalization(a,d,energ))).*energ
    return ρ
end

# A manual search for energy minima with plot (admissability check ignored)
using Plots
Plots.plot(x->enplothelp(x),0.8371,0.8374,color=:black,label="E(R)",xlabel="R",ylabel="E(R)",grid=false)
# Use Optim.jl or other comparable packages to compute the minimizer to desired accuracy and check admissibility
Rmin = 0.8372415
Plots.scatter!((Rmin,enplothelp(Rmin)),label="computed minimizer")
enplothelp(Rmin)
# Obtain the equilibrium measure
ρcfs = computemeas(Rmin,op1,op2,d,α,β,a,b,M)
ρfunc(x) = (1-x^2)^a*(Jacobi(a,b)[QuadraticMap(),1:length(ρcfs)]*ρcfs)[x]

# Plot obtained measure (radial cut from r = (0,R))
plot(x->(1/(Rmin^d)*ρfunc(x/Rmin)),0:0.001:Rmin,grid=false,thickness_scaling=1.2,xlabel="R",ylabel="ρ",color=:black,label="computed ρ")
