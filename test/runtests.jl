using EquilibriumMeasures, IntervalSets
# using Plots

μ = equilibriummeasure(x -> x^2)
a,b = first(axes(μ,1)),last(axes(μ,1))
xx = range(a,b; length=20)
view(μ,xx)

plot(xx, getindex.(Ref(μ), xx))

@time μ[0.1]

μ.args[1][0.1,
μ.args[2].datasize

μ[xx]

OrthogonalPolynomialsQuasi, ContinuumArrays, IntervalSets, FillArrays, ArrayLayouts, LazyArrays
using ForwardDiff
Base.floatmin(::Type{<:ForwardDiff.Dual}) = floatmin(Float64)

a = ForwardDiff.Dual(sqrt(2),1)
a = 2.0

import LazyArrays: paddeddata






using ForwardDiff
mom(V, ForwardDiff.Dual(2.0,1))



function mom(V,a)
    x = Inclusion(-a..a)
    T = ChebyshevT{typeof(a)}()[x/a,:]
    U = ChebyshevU{typeof(a)}()[x/a,:]
    cfs = ldiv(T , V.(x))
    D = apply(*, Derivative(axes(T,1)), T).args[2]
    UT = ldiv(U, T)
    Vp = UT \ (D * cfs)
    
end


Ta = (V,a) -> begin
    w, U = ChebyshevUWeight{typeof(a)}(), ChebyshevU{typeof(a)}()
    x = Inclusion(-a..a)
    H = inv.(x .- x')
    # basis on -a..a
    T = ChebyshevT{typeof(a)}()[x/a,:]
    wU = (w .* U)[x/a,:] 
    # Hilbert
    App =  applied(*, H, wU)
    HwU = materialize(App)
    LT = ldiv(T, HwU)
    Hd = LT[2:end,:]
    Vf = T* ldiv(T,V.(x))
    Vp = diff(Vf)
    w = wU * ldiv(Hd, ldiv(T, Vp)[2:end])
end



a = 2.0; x = Inclusion(-a..a)
@enter view(ChebyshevT{typeof(a)}(),x/a,:)
Ta(V, a)
using DualNumbers


a =  dual(1.0,1)


Ta = (V,a) -> begin
    w, U = ChebyshevUWeight{typeof(a)}(), ChebyshevU{typeof(a)}()
    x = Inclusion(-a..a)
    H = inv.(x .- x')
    # basis on -a..a
    T = ChebyshevT{typeof(a)}()[x/a,:]
    wU = (w .* U)[x/a,:] 
    # Hilbert
    App =  applied(*, H, wU)
    HwU = materialize(App)
    LT = ldiv(T, HwU)
    Hd = LT[2:end,:]
    Vf = T* ldiv(T,V.(x))
    Vp = diff(Vf)
    w = wU * ldiv(Hd, ldiv(T, Vp)[2:end])
end

using QuasiArrays
import ContinuumArrays: demap, ApplyQuasiArray
A,B = T, HwU; A,B = demap(A),demap(first(B.args)); B = ApplyQuasiArray(B);
@which apply(\,A,B)
ldiv(T, HwU)





V = x->x^2; a = 2.0

@code_warntype Ta(x->x^2, 2.0)

)
import OrthogonalPolynomialsQuasi: adaptivetransform_ldiv
V = x -> x^2; T = Chebyshev(); @which adaptivetransform_ldiv(T,V.(axes(T,1)))
x = Inc
x
ForwardDiff.derivative(a -> begin
    T = Chebyshev{typeof(a)}()
    x = axes(T,1); T \ V.(x) 
end, 2.0)

T = Chebyshev(); x = axes(T,1)
V = x -> x^2; Vx = V.(x)

@which T \ Vx

@time ForwardDiff.partials(sum(Ta(V, ForwardDiff.Dual(2.0,1)))/2)[1]
Sx = x -> sum(Ta(V,x))/2
@time ForwardDiff.derivative(Sx, 2.0)
sum(w)/2

equilibriummeasure(x -> log.(abs.(x .- x')))
equilibriummeasure(x -> abs.(x .- x').^α)

using ForwardDiff


T = Chebyshev()
U = Ultraspherical(1) 
x = axes(T,1)
V = x.^2
H = inv.(x .- x')
wU = UltrasphericalWeight(1) .* Ultraspherical(1)



sum(w)







Ta = a -> begin
    x = Inclusion(-a..a)
    H = inv.(x .- x')
    T = Chebyshev{typeof(a)}()[x/a,:]
    Tn = Chebyshev{typeof(a)}()[x/a,2:5]
    wU = (w .* U)[x/a,:]
    (T \ (H*wU))[2,1] \ (Tn \ x)[1]
end

N = a -> Ta(a)[1] - sqrt(2)/π


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
