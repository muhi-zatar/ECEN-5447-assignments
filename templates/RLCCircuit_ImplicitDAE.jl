cd(@__DIR__)
using Pkg
Pkg.activate(".")
#Pkg.resolve() # make sure Manifest matches Project
#Pkg.instantiate() # install missing dependencies

# Import packages
using DifferentialEquations
using Plots
using LinearAlgebra
using Sundials

# Define function to be used by integration method  
function implicitDAE(out, du, u, p, t)
    # Rename locally for readability
    R, L, C = p
    vc, il, vg, ig, vr, ir, vl, ic = u

    # Define the system of equations
    # NOTE 1: out[] indexing is different than u,du indexing
    # NOTE 2: Solver will try to make every LHS zero
    out[1] = -du[1] + (1/C)*ic
    out[2] = -du[2] + (1/L)*vl
    out[3] = -vg + 1
    out[4] = -ig - il
    out[5] = -vr + vc
    out[6] = -ir + (1/R)*vr
    out[7] = -vl + vg - vc
    out[8] = -ic + il - ir
end

# Define initial conditions
u0 = [0.8, 5.0, 1.0, -5.0, 0.8, 0.16, 0.2, 4.84]
du0 = [48.4, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Define length of simulation
tspan = (0.0, 5.0)

# Define parameters
R = 5.0; L = 0.1; C = 0.1;
p = [R, L, C]

# Define which variables have differential terms 
diff_vars = [true, true, false, false, false, false, false, false]

# Build problem
prob = DAEProblem(implicitDAE, du0, u0, tspan, p; 
                  differential_vars=diff_vars)
# Run simulation
sol = solve(prob, IDA())

# Plot
plot(sol, idxs=(0,[1,2]), # (vc, il) vs time
    title="RLC Circuit - Implicit DAE (Residual)", 
    labels=["vc" "il"]
    )


# ------------------------------ EXTRAS ------------------------------
# Check IC consistency 
out = Vector{Float64}(undef,8); # return object for implicitDAE()
implicitDAE(out, du0, u0, p, 0); # evaluate function once
out # look at residual vector

# Solve system for one time step using an initializer to find consistent ICs
u0_guess = [0.8, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
du0_guess = zeros(8)
tspan_init = (0.0, 0.0)
prob_init = DAEProblem(
    implicitDAE, du0_guess, u0_guess, tspan_init, p; 
    differential_vars=diff_vars, 
    initializealg=BrownFullBasicInit()
    )
sol_init = solve(prob_init, DFBDF())
sol_init.u[1]

# Before solving with initializer, put this inside the function to find du
# (For some reason, it is hard to access du0 after the solve)
if t==0
    print("### [time $t] du:")
    println(du)
    print("### [time $t] u:")
    println(u)
    print("### [time $t] out:")
    println(out)
end

