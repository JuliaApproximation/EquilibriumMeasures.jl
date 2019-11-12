using OrthogonalPolynomialsQuasi, ContinuumArrays, IntervalSets

x = Inclusion(-sqrt(2)..sqrt(2))
H = inv.(x .- x')

w, U = UltrasphericalWeight(1), Ultraspherical(1)
wU = (w .* U)[x/sqrt(2),:]
T = Chebyshev()[x/sqrt(2),:]

T \ (H*wU)
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
