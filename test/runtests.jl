using OrthogonalPolynomialsQuasi, ContinuumArrays, IntervalSets, FillArrays

x = Inclusion(-sqrt(2)..sqrt(2))


equilibriummeasure(x -> log.(abs.(x .- x')))
equilibriummeasure(x -> abs.(x .- x').^α)

using ForwardDiff

w, U = UltrasphericalWeight(1), Ultraspherical(1)

Ta = a -> begin
    x = Inclusion(-a..a)
    H = inv.(x .- x')
    T = Chebyshev{typeof(a)}()[x/a,:]
    Tn = Chebyshev{typeof(a)}()[x/a,2:5]
    wU = (w .* U)[x/a,:]
    ((T \ (H*wU))[2:end,:] \ [Tn \ x; Zeros{typeof(a)}(∞)])
end

Ta = a -> begin
    x = Inclusion(-a..a)
    H = inv.(x .- x')
    T = Chebyshev{typeof(a)}()[x/a,:]
    Tn = Chebyshev{typeof(a)}()[x/a,2:5]
    wU = (w .* U)[x/a,:]
    (T \ (H*wU))[2,1] \ (Tn \ x)[1]
end

N = a -> Ta(a)[1] - sqrt(2)/π
Base.floatmin(::Type{<:ForwardDiff.Dual}) = floatmin(Float64)

ForwardDiff.derivative(N,1.4)

N(1.4+0.00001)-N(1.4)


Tn = Chebyshev()[x/sqrt(2),2:5]



using Plots
g = range(-1.4,1.4; length = 100)
plot(g, Ta(1.4)[g])






Vp = T * [Tn \ (2x); zeros(∞)]





(Tn * (Tn \ x))[0.1]


T = Chebyshev()
x = axes(T,1)
T[:,1:2] \ x



# need to add transforms to OrthogonalPolynomialsQuasi
g = range(-sqrt(2),sqrt(2); length=3)
Vf = x -> x^2
V = T * [T[g,1:3] \ Vf.(g); zeros(∞)]
D = Derivative(x)
(U[x/sqrt(2),:] \ T) \  (U[x/sqrt(2),:] \ (D*V))

grid(Chebyshev())

x.^2
findfirst(isequal(0.1),2x)
2x

x = axes(U,1)


(H*wU)


using EquilibriumMeasures, ApproxFun, Plots




μ = equilibriummeasure(x -> x^2)
equilibriummeasure(x -> x, 0..1; bounded=:right)
plot(μ)
