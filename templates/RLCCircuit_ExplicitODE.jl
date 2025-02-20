cd(@__DIR__)
using Pkg
Pkg.activate(".")
#Pkg.resolve() # make sure Manifest matches Project
#Pkg.instantiate() # install missing dependencies

# Import packages
using DifferentialEquations
using Plots

# Define function to be used by integration method 
function explicitODE(du, u, p, t)
    # Rename locally for readability
    R, L, C = p
    vc, il = u

    # Define system of equations
    du[1] = (1/C) * (il - (1/R)*vc) # d(vc)/dt 
    du[2] = (1/L) * (1 - vc)        # d(il)/dt 
end

# Define initial conditions
u0 = [0.8, 5.0] # [vc, il]

# Define length of simulation
tspan = (0.0, 5.0)

# Define parameters
R = 5.0
L = 0.1
C = 0.1
p = [R, L, C]

# Build problem
prob = ODEProblem(explicitODE, u0, tspan, p)

# Run simulation
sol = solve(prob, RK4())

# Plot
plot(sol, title = "RLC Circuit - Explicit ODE", labels = ["vc" "il"])


# ------------------------------ EXTRAS ------------------------------
# Plot the phase portrait of the two state variables
plot(sol, 
    idxs=(1,2),  # 0 = time, 1 = 1st state, 2 = 2nd state, ...
    title="Phase Portrait", 
    xlabel="Capacitor Voltage", 
    ylabel="Inductor Current",
    label=:none,
    linewidth=3
    )
scatter!([u0[1]], [u0[2]], label="Initial Condition")
