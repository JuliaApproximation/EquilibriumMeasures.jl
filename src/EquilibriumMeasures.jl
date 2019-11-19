
module EquilibriumMeasures
using Base, OrthogonalPolynomialsQuasi, ContinuumArrays, ForwardDiff

import ForwardDiff: derivative

export equilibriummeasure

function _equilibriummeasure(V, a)
    x = Inclusion(-a..a)
    T = ChebyshevT{typeof(a)}()[x/a,:]
    wU = (ChebyshevUWeight{typeof(a)}() .* ChebyshevU{typeof(a)}())[x/a,:]
    U = ChebyshevU{typeof(a)}()[x/a,:]
    # Operators
    H = T \ (inv.(x .- x') * wU) # Hilbert wU -> T
    D = U \ (Derivative(axes(T,1)) * T) # Derivative T -> U
    UT = U \ T  # Converion T -> U
    V_cfs = T \ V.(x)
    Vp_cfs = UT \ (D * V_cfs)
    wU * (H[2:end,:] \ Vp_cfs[2:end])
end

# Intentionally hide type for compile time
struct EquilibriumMeasureMoment
    V
end

(E::EquilibriumMeasureMoment)(a) = sum(_equilibriummeasure(E.V,a))/2 - 1


function equilibriummeasure(V; a = 2.0, maxiterations=1000)
    μ = EquilibriumMeasureMoment(V)
    for k=1:maxiterations    
        an = a - derivative(μ, a)\μ(a)
        an ≈ a && return _equilibriummeasure(V,a)
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
