
module EquilibriumMeasures
using Base, OrthogonalPolynomialsQuasi


# using Base, Compat, DualNumbers, ApproxFun, SingularIntegralEquations

# export equilibriummeasure,hilbertinverse




# ## hilbertinverse
# # returns a function whose hilbert transform is the given function

# function hilbertinverse(u::Fun;tolerance=100eps(),bounded=:none)
# cfs=coefficients(u,Chebyshev)
# if abs(first(cfs)) < tolerance
#     # no singularity
#     # invert Corollary 5.6 of Olver&Trogdon
#     cfs=-coefficients(cfs[2:end], Ultraspherical(1,domain(u)),Chebyshev(domain(u)))
#     Fun(JacobiWeight(.5,.5,Chebyshev(domain(u))),cfs)
# else
#     # invert Corollary 5.10 of Olver&Trogdon
#     if bounded==:left
#         cfs=[cfs[1];cfs[1];.5*cfs[2:end]]
#     elseif bounded==:right
#         cfs=[-cfs[1];cfs[1];.5*cfs[2:end]]
#     else
#         cfs=[0.;cfs[1];.5*cfs[2:end]]
#     end
#     cfs=coefficients(cfs,ChebyshevDirichlet{1,1}(domain(u)),Chebyshev(domain(u)))
#     Fun(JacobiWeight(-.5,-.5,Chebyshev(domain(u))),cfs)
# end
# end



# # this must be zero for the corrects support

# function emmoment(f::Fun{<:Chebyshev};bounded=:none)
# cfs=pad(f.coefficients,2)
# d=domain(f)
# a,b = endpoints(d)
# if bounded==:left
#     (a-b)/4*cfs[1]
# elseif bounded==:right
#     (b-a)/4*cfs[1]
# else
#     (b-a)/8*cfs[2]
# end
# end

# function F(V,a,b)
# f=Fun(V,a..b)'

# [f.coefficients[1],emmoment(f)-1]
# end


# equilibriummeasuresupport(V,ab::Vector;opts...) =
# equilibriummeasuresupport(V,convert(Domain,ab);opts...)
# function equilibriummeasuresupport(V,ab=(-1.0..1.0);maxiterations=100,bounded=:none)
# if bounded == :left
#     return equilibriummeasuresupportbounded(false,V,ab;maxiterations=maxiterations)
# elseif bounded == :right
#     return equilibriummeasuresupportbounded(true,V,ab;maxiterations=maxiterations)
# end

# a,b = endpoints(ab)

# # Newton iteration, using dual numbers
# for k=1:maxiterations
#     F1 = F(V,dual(a,1.),b)
#     p = realpart.(F1)
#     J=[epsilon.(F1) epsilon.(F(V,a,dual(b,1)))]

#     an,bn=sort([a,b]-J\p)
#     if isapprox(an,a) && isapprox(bn,b)
#         return (an..bn)
#     else
#         a,b=an,bn
#     end
# end

# warn("maxiterations reached")
# a..b
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
