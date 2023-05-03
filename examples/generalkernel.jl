##
# consider instead of the Log kernel we want the equilibrium measure with kernel H[sin(x-t)]
##

using EquilibriumMeasures, MultivariateOrthogonalPolynomials, ClassicalOrthogonalPolynomials, FillArrays

K = (x,y) -> sin(x-t)

RectPolynomial(Fill(Legendre(), 2))

expand(