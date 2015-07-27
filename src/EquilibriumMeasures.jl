module EquilibriumMeasures
    using Base, Compat, DualNumbers, ApproxFun, SingularIntegralEquations

export equilibriummeasure

# We need to implement some functionality for the ApproxFun constructor to work with dual numbers
Base.sinpi(x::Dual)=sin(π*x)
ApproxFun.real{T}(::Type{Dual{T}})=Dual{ApproxFun.real(T)}
ApproxFun.plan_chebyshevtransform{D<:Dual}(v::Vector{D})=ApproxFun.plan_chebyshevtransform(real(v))
ApproxFun.chebyshevtransform{D<:Dual}(v::Vector{D},plan...)=dual(chebyshevtransform(real(v),plan...),chebyshevtransform(epsilon(v),plan...))
ApproxFun.chop!(f::Fun,d::Dual)=chop!(f,real(d))




## hilbertinverse
# returns a function whose hilbert transform is the given function

function hilbertinverse(u::Fun;tolerance=100eps(),bounded=:none)
    cfs=coefficients(u,Chebyshev)
    if abs(first(cfs)) < tolerance
        # no singularity
        # invert Corollary 5.6 of Olver&Trogdon
        cfs=-coefficients(cfs[2:end],Ultraspherical{1},Chebyshev)
        Fun(cfs,JacobiWeight(.5,.5,Chebyshev(domain(u))))
    else
        # invert Corollary 5.10 of Olver&Trogdon
        if bounded==:left
            cfs=[cfs[1],cfs[1],.5*cfs[2:end]]
        elseif bounded==:right
            cfs=[-cfs[1],cfs[1],.5*cfs[2:end]]
        else
            cfs=[0.,cfs[1],.5*cfs[2:end]]
        end
        cfs=coefficients(cfs,ChebyshevDirichlet{1,1},Chebyshev)
        Fun(cfs,JacobiWeight(-.5,-.5,Chebyshev(domain(u))))
    end
end



# this must be zero for the corrects support

function emmoment(f::Fun{Chebyshev};bounded=:none)
    cfs=pad(f.coefficients,2)
    if bounded==:left
        (a-b)/4*cfs[1]
    elseif bounded==:right
        (b-a)/4*cfs[1]
    else
        (b-a)/8*cfs[2]
    end
end

function F(V,a,b)
    f=Fun(V,Interval(a,b))'

    [f.coefficients[1],emmoment(f)-1]
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
