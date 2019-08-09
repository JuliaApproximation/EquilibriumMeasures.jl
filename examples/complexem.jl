using ApproxFun, DualNumbers

import EquilibriumMeasures: emmoment

Vp = x -> -im*x^2


# Vp = x -> im
a,b = -1.1+0.5im,1+im
a,b = 1,im
N = (a,b) -> (f = Fun(Vp, Segment(a,b)); [f.coefficients[1],emmoment(f)-1])

J=[epsilon.(N(dual(a,1),b)) epsilon.(N(a,dual(b,1)))]
a,b = [a,b]-J\N(a,b)

plot(Segment(a,b))