
module EquilibriumMeasures
using Base, OrthogonalPolynomialsQuasi, ContinuumArrays, ForwardDiff, IntervalSets, DomainSets, StaticArrays

import ForwardDiff: derivative, gradient, jacobian

export equilibriummeasure

Base.floatmin(::Type{<:ForwardDiff.Dual}) = floatmin(Float64)

function _equilibriummeasure(V, a, b)
    Typ = float(promote_type(typeof(a), typeof(b)))
    T = ChebyshevT{Typ}()
    x = Inclusion(a..b)
    y = affine(x, axes(T,1)) # affine map from a..b to -1..1
    T̃ = T[y,:]
    wŨ = (ChebyshevUWeight{Typ}() .* ChebyshevU{Typ}())[y,:]
    Ũ = ChebyshevU{Typ}()[y,:]
    # Operators
    H = T̃ \ (inv.(x .- x') * wŨ) # Hilbert wU\-> T̃
    D = Ũ \ (Derivative(axes(T̃,1)) * T̃) # Derivative T̃ -> Ũ
    ŨT̃ = Ũ \ T̃  # Converion T̃ -> Ũ
    V_cfs = T̃ \ V.(x)
    Vp_cfs = ŨT̃ \ (D * V_cfs)
    Vp_cfs[1], wŨ * (H[2:end,:] \ Vp_cfs[2:end])
end


# Intentionally hide type for compile time
struct EquilibriumMeasureMoment
    V
end

function (E::EquilibriumMeasureMoment)(a) 
    c0,μ = _equilibriummeasure(E.V, a...)
    SVector(c0, sum(μ)/2 - 1)
end


function equilibriummeasure(V; a = SVector(-1.0,1.0), maxiterations=1000)
    μ = EquilibriumMeasureMoment(V)
    for k=1:maxiterations    
        an = a - jacobian(μ, a)\μ(a)
        an ≈ a && return _equilibriummeasure(V, a...)[2]
        a = an
    end
    error("Max its")
end



# function equilibriummeasuresupport(V,ab=(-1.0..1.0);maxiterations=100,bounded=:none)
#     if bounded == :left
#         return equilibriummeasuresupportbounded(false,V,ab;maxiterations=maxiterations)
#     elseif bounded == :right
#         return equilibriummeasuresupportbounded(true,V,ab;maxiterations=maxiterations)
#     end

#     a,b = endpoints(ab)

#     # Newton iteration, using dual numbers
#     for k=1:maxiterations
#         F1 = F(V,dual(a,1.),b)
#         p = realpart.(F1)
#         J=[epsilon.(F1) epsilon.(F(V,a,dual(b,1)))]

#         an,bn=sort([a,b]-J\p)
#         if isapprox(an,a) && isapprox(bn,b)
#             return (an..bn)
#         else
#             a,b=an,bn
#         end
#     end

#     warn("maxiterations reached")
#     a..b
# end




# Fbounded(s::Bool,V,a,b) = emmoment(Fun(V,(a..b))'; bounded= s ? (:right) : (:left))-1

# function equilibriummeasuresupportbounded(s::Bool,V,ab=(-1.0..1.0);maxiterations=100)
# a,b=endpoints(ab)
# for k=1:maxiterations
#     F1=s ? Fbounded(s,V,a,dual(b,1.)) : Fbounded(s,V,dual(a,1.),b)

#     p=realpart(F1)
#     J=epsilon(F1)

#     if s
#         bn=b-J\p
#         if isapprox(bn,b)
#             return (a..bn)
#         else
#             b=bn
#         end
#     else
#         an=a-J\p
#         if isapprox(an,a)
#             return (an..b)
#         else
#             a=an
#         end
#     end
# end
# warn("maxiterations reached")
# a..b
# end

# # we assume that the zeroth coefficient is zero
# equilibriummeasure(V,gs...;bounded=:none,opts...) =
# -hilbertinverse(Fun(V,equilibriummeasuresupport(V,gs...;bounded=bounded,opts...))';bounded=bounded)/(2π)

end #module
