module EquilibriumMeasures
    using Base, Compat, DualNumbers, ApproxFun, SingularIntegralEquations

export equilibriummeasure

# We need to implement some functionality for the ApproxFun constructor to work with dual numbers
Base.sinpi(x::Dual)=sin(π*x)
ApproxFun.real{T}(::Type{Dual{T}})=Dual{ApproxFun.real(T)}
ApproxFun.plan_chebyshevtransform{D<:Dual}(v::Vector{D})=ApproxFun.plan_chebyshevtransform(real(v))
ApproxFun.chebyshevtransform{D<:Dual}(v::Vector{D},plan...)=dual(chebyshevtransform(real(v),plan...),chebyshevtransform(epsilon(v),plan...))
ApproxFun.chop!(f::Fun,d::Dual)=chop!(f,real(d))


# this must be zero for the corrects support
function F(V,a,b)
    f=Fun(V,Interval(a,b))'
    [f.coefficients[1],(b-a)/8*f.coefficients[2]-1]
end



function equilibriummeasuresupport(V,ab=Interval();maxiterations=100)
    a,b=sort([ab.a,ab.b])

    # Newton iteration, using dual numbers
    for k=1:maxiterations
        F1=F(V,dual(a,1.),b)
        p=real(F1)
        J=[epsilon(F1) epsilon(F(V,a,dual(b,1)))]

        an,bn=sort([a,b]-J\p)
        if isapprox(an,a) && isapprox(bn,b)
            return Interval(an,bn)
        else
            a,b=an,bn
        end
    end

    warn("maxiterations reached")
    Interval(a,b)
end
# we assume that the zeroth coefficient is zero
equilibriummeasure(V,gs...;opts...)=-hilbertinverse(Fun(V,equilibriummeasuresupport(V,gs...;opts...))';tolerance=0.1)/(2π)

end #module
