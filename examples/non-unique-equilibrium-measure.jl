using EquilibriumMeasures, StaticArrays
using Plots

"""
This script finds 9 solutions for the equilibrium measure of which 3 are feasible. 
"""

# Potential
V = x -> (x-3)*(x-2)*(1+x)*(2+x)*(3+x)*(2x-1)/20

# Initial guess
ic = SVector(-1,1)

# First solution
a = equilibriummeasure(V; a=ic)

# Keep first solution in an array
solns = [a]

for no_sols = 1:2
    
    # Compute next solution with the same initial guess but deflating
    # all other known solutions
    b = equilibriummeasure(V; a=ic, knownsolutions=solns, dampening=0.5)

    # Append array of known solutions with new solution
    push!(solns, b)

end

# Try different damping
b = equilibriummeasure(V; a=ic, knownsolutions=solns, dampening=0.3)
push!(solns, b)

# Try different initial guesses
ic = SVector(2,3)
b = equilibriummeasure(V; a=ic, knownsolutions=solns, dampening=0.5)
push!(solns, b)

ic = SVector(1,3)
b = equilibriummeasure(V; a=ic, knownsolutions=solns, dampening=0.5)
push!(solns, b)

ic = SVector(-1,4)
b = equilibriummeasure(V; a=ic, knownsolutions=solns, dampening=0.5)
push!(solns, b)

b = equilibriummeasure(V; a=ic, knownsolutions=solns, dampening=0.4)
push!(solns, b)

b = equilibriummeasure(V; a=ic, knownsolutions=solns, dampening=0.5)
push!(solns, b)

# Plot solutions
p = plot(_equilibriummeasure(V, solns[1]...)[2])
if length(solns) > 1
    for i = 2:length(solns)
        p = plot!(_equilibriummeasure(V, solns[i]...)[2], legend=:top)
    end
end

display(p)
# savefig(p, "em.pdf")