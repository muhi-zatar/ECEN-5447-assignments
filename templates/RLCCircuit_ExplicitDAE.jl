cd(@__DIR__)
using Pkg
Pkg.activate(".")
#Pkg.resolve() # make sure Manifest matches Project
#Pkg.instantiate() # install missing dependencies

# Import packages
using DifferentialEquations
using Plots
using LinearAlgebra

# Define function to be used by integration method  
function explicitDAE(du, u, p, t)
    # Rename locally for readability
    R, L, C = p
    vc, il, vg, ig, vr, ir, vl, ic  = u

    # Define the system of equations
    du[1] = (1/C)*ic        # d/dt (vc) != 0
    du[2] = (1/L)*vl        # d/dt (il) != 0
    du[3] = -vg + 1         # d/dt (vg) = 0
    du[4] = -ig - il        # d/dt (ig) = 0
    du[5] = -vr + vc        # d/dt (vr) = 0
    du[6] = -ir + (1/R)*vr  # d/dt (ir) = 0
    du[7] = -vl + vg - vc   # d/dt (vl) = 0
    du[8] = -ic + il - ir   # d/dt (ic) = 0
end

# Build mass matrix
M = Diagonal([ones(2); zeros(6)])

# Build function 
explicitDAE_M = ODEFunction(explicitDAE, mass_matrix = M)

# Define initial conditions
u0 = [0.8, 5.0, 1.0, -5.0, 0.8, 0.16, 0.2, 4.84]

# Define length of simulation
tspan = (0.0, 5.0)

# Define parameters
R = 5.0
L = 0.1
C = 0.1
p = [R, L, C]

# Build problem
prob = ODEProblem(explicitDAE_M, u0, tspan, p)

# Run simulation
sol = solve(prob, Rodas5());

# Plot
plot(sol, idxs=(0,[1,2]), # 1st & 2nd states (vc, il) vs. time
    title="RLC Circuit - Explicit DAE (Mass Matrix)", 
    labels=["vc" "il"]
    )


# ------------------------------ EXTRAS ------------------------------
# Check IC consistency 
du0 = Vector{Float64}(undef,8); # return object for explicitDAE()
explicitDAE(du0, u0, p, 0); # evaluate function once
du0 # look at derivative vector 