module EquilibriumMeasures
using Base, ClassicalOrthogonalPolynomials, ContinuumArrays, ForwardDiff, IntervalSets, DomainSets, StaticArrays,
    FillArrays, BlockArrays, LazyBandedMatrices, FastTransformsForwardDiff, SingularIntegrals

import ForwardDiff: derivative, gradient, jacobian

import LinearAlgebra: dot
import IntervalSets: mean

export equilibriummeasure, _equilibriummeasure, deflation_inner_products, deflation_deriv, deflation_op, deflation_scale

Base.floatmin(::Type{<:ForwardDiff.Dual}) = floatmin(Float64)

function _equilibriummeasure(V, ab::SVector{2})
    a,b = ab
    Typ = float(eltype(ab))
    T = ChebyshevT{Typ}()
    x = Inclusion(a..b)
    y = affine(x, axes(T,1)) # affine map from a..b to -1..1
    T̃ = T[y,:]
    wŨ = Weighted(ChebyshevU{Typ}())[y,:]
    Ũ = ChebyshevU{Typ}()[y,:]
    # Operators
    H = T̃ \ (inv.(x .- x') * wŨ) # Hilbert wU\-> T̃
    D = Ũ \ (Derivative(axes(T̃,1)) * T̃) # Derivative T̃ -> Ũ
    ŨT̃ = Ũ \ T̃  # Converion T̃ -> Ũ
    V_cfs = T̃ \ V.(x)
    Vp_cfs = ŨT̃ \ (D * V_cfs)
    Vp_cfs[1], wŨ * (H[2:end,:] \ Vp_cfs[2:end])/2
end

function _equilibriummeasure(V, ad::SVector{4})
    a,b,c,d = ad
    Typ = float(eltype(ad))
    T1,T2 = chebyshevt(a..b),chebyshevt(c..d)
    U1,U2 = chebyshevu(a..b),chebyshevu(c..d)
    T = PiecewiseInterlace(T1, T2)
    U = PiecewiseInterlace(U1, U2)
    W = PiecewiseInterlace(Weighted(U1), Weighted(U2))
    x = axes(W,1)
    # Operators
    D = U \ (Derivative(x) * T) # Derivative T -> U
    UT = U \ T  # Converion T -> U
    V_cfs = T \ V.(x)
    Vp_cfs = UT \ (D * V_cfs)

    H = T \ (inv.(x .- x') * W) # Hilbert W -> T

    # add extra columns for constants
    H̃ = BlockHcat(Eye((axes(H,1),))[:,Block(1)], 2H)
    c = H̃ \ Vp_cfs
    c[Block(1)], W * c[Block.(2:∞)]
end


# Intentionally hide type for compile time
struct DFunction <: Function
    f
end
(f::DFunction)(x::T) where T = f.f(x)::T
Base.convert(::Type{DFunction}, f::DFunction) = f
Base.convert(::Type{DFunction}, f::Function) = DFunction(f)


struct EquilibriumMeasureMoment
    V::DFunction
end

_logterms(V, μ) = ()
_logterms(V, μ, d) = ()
function _logterms(V, μ, d1, d2)
    x = axes(μ,1)
    z1,z2 = mean(d1),mean(d2)
    (2*(log.(abs.(z1 .- x'))*μ) - V(z1) - 2*(log.(abs.(z2 .- x'))*μ) + V(z2),)
end

function (E::EquilibriumMeasureMoment)(a)
    c0,μ = _equilibriummeasure(E.V, a)
    SVector(c0..., sum(μ) - 1, _logterms(E.V, μ, components(axes(μ,1).domain)...)...)
end

equilibriummeasure(V; a = SVector(-1.0,1.0), maxiterations=1000, knownsolutions=[], power=2, shift=1.0, dampening=1.0, returnendpoint=false) =
    deflate_equilibriummeasure(EquilibriumMeasureMoment(DFunction(V)), a, maxiterations, knownsolutions, power, shift, dampening, returnendpoint)
function deflate_equilibriummeasure(μ, a, maxiterations, knownsolutions, power, shift, dampening, returnendpoint)
    num_found_sols = length(knownsolutions)
    for k=1:maxiterations
        update = - jacobian(μ, a)\μ(a)
        if num_found_sols > 0
            update = deflation_scale(update, a, knownsolutions, power, shift) * update
        end
        an = a + dampening*update
        if an ≈ a
            # improve accuracy a bit more
            for _=1:3
                update = - jacobian(μ, a)\μ(a)
                if num_found_sols > 0
                    update = deflation_scale(update, a, knownsolutions, power, shift) * update
                end
                a = a + dampening*update
            end
            if returnendpoint
                return  _equilibriummeasure(μ.V, a)[2], a
            else
                return  _equilibriummeasure(μ.V, a)[2]
            end
        end
        a = an
    end
    error("Max its")
end

# Deflation code

function  deflation_inner_products(state::AbstractVector, sol)
    num_sols = length(sol)
    T = promote_type(eltype(state), eltype(eltype(sol)))

    inner_products = zeros(T, num_sols)
    for iter = 1:num_sols
        inner_products[iter] = dot(state.-sol[iter], state.-sol[iter])
    end
    inner_products
end

function deflation_op(state::AbstractVector, sol, power::Real, shift::Real)
    num_sols = length(sol)
    T = promote_type(eltype(state), eltype(eltype(sol)))
    power = T(power); shift = T(shift)

    m = one(T)
    inner_products = deflation_inner_products(state, sol)

    for iter = 1:num_sols
        normsq = inner_products[iter]
        factor = one(T)/normsq^(power/2) + shift
        m = m*factor
    end
    m, inner_products
end

function deflation_deriv(update::AbstractVector, state::AbstractVector, sol::Vector, power::Real, shift::Real)
    m, inner_products = deflation_op(state, sol, power, shift)
    num_sols = length(sol)
    T = promote_type(eltype(state), eltype(eltype(sol)))
    power = T(power); shift = T(shift)

    dm = zero(T)
    state_dot_update = dot(state, update)
    for iter = 1:num_sols
        scale = m/deflation_op(state, sol[iter], power, shift)[1]
        deriv = -power*state_dot_update/inner_products[iter]^((power + 1)/2)
        dm += scale*deriv
    end
    m, dm
end

function deflation_scale(update::AbstractVector, state::AbstractVector, sol::Vector, power::Real, shift::Real)
    T = promote_type(eltype(state), eltype(eltype(sol)))
    m, dm = deflation_deriv(update, state, sol, power, shift)
    minv = one(T)/m
    tau = one(T) + minv*dm/(one(T) - minv*dm)
    tau
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
