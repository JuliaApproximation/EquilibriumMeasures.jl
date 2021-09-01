# EquilbriumMeasures
Calculate equilibrium measures from potentials

[![Build Status](https://github.com/JuliaApproximation/EquilibriumMeasures.jl/workflows/CI/badge.svg)](https://github.com/JuliaApproximation/EquilibriumMeasures.jl/actions)
[![codecov](https://codecov.io/gh/JuliaApproximation/EquilibriumMeasures.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/EquilibriumMeasures.jl)


In potential theory, the equilibrium measure is the limiting distribution of charges in a potential field that minimises energy, see [Saff and Totic 1997](https://www.springer.com/gp/book/9783540570783) and [Deift 1999](https://bookstore.ams.org/cln-3). These distributions also describe the distribution of
1. Zeros of orthogonal polynomials
2. Fekete (near optimal interpolation) points
3. Eigenvalues of random matrices
This package facilitates their computation, and is based on [RHPackage](https://github.com/dlfivefifty/RHPackage).

```julia
using EquilibriumMeasures, Plots
plot(equilibriummeasure(x -> x^2); label="Semicircle", legend=:topleft)
plot!(equilibriummeasure(x -> x^4); label="x^4")
plot!(equilibriummeasure(x -> x^100); label="x^100")
plot!(equilibriummeasure(x -> exp(x) - x); label="exp(x) - x")
```
<img src=https://github.com/JuliaApproximation/EquilibriumMeasures.jl/raw/master/images/ems.svg>