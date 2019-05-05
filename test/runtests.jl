using EquilibriumMeasures, ApproxFun, Plots

μ = equilibriummeasure(x -> x^2)
equilibriummeasure(x -> x, 0..1; bounded=:right)
plot(μ)
