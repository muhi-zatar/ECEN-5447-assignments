cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.resolve()
Pkg.instantiate()

# Import packages
using DifferentialEquations
using LinearAlgebra
using Plots
using PowerSystems
using PowerFlows
const PSY = PowerSystems
const PF = PowerFlows

# -----------------------------------------------------------------------------------------
# Section 1: Build system and run powerflow 
# -----------------------------------------------------------------------------------------

# Load raw file into Sienna
dir_raw = "data/ThreeBusMultiLoad.raw"
sys = PSY.System(dir_raw)

# Run power flow
pf_result = PF.solve_powerflow(ACPowerFlow{}(), sys)
PSY.set_units_base_system!(sys, "SYSTEM_BASE")

# Get power flow output: complex voltage of each bus terminal 
# (Sienna exports these in pu and rad)
v = pf_result["bus_results"].Vm # [pu-V]
θ = pf_result["bus_results"].θ  # [rad]

# Get power flow output: real (P) and reactive (Q) power injected at each bus
# (Sienna exports these in MW/MVA, so adjust by system base)
P = pf_result["bus_results"].P_net / PSY.get_base_power(sys) # [pu(MW)]
Q = pf_result["bus_results"].Q_net / PSY.get_base_power(sys) # [pu(MVar)]


# -----------------------------------------------------------------------------------------
# Section 2: Calculate relevant quantities needed for running a time domain simulation
# -----------------------------------------------------------------------------------------

##### TODO: Using v,θ,P,Q to find the parameters and initial conditions

# To get you started, here are the branch parameters
PSY.show_components(Line, sys, [:r, :x, :b, :arc])

# Find complex power
S = P + im * Q

# Find complex voltage
V = v .* (cos.(θ) .+ im .* sin.(θ))

# Find complex current
I = conj(S ./ V)

# DQ components are the real and imaginary components of the complex expressions
# Node voltages are differential variables, so power flow solution provides initial condition
vd_init = real(V)
vq_init = imag(V)

# Node currents are parameters
id = real(I)
iq = imag(I)

# Define some arbitrary generator reactances
x1 = 1.0
x2 = 1.0

# Find load impedance
Z = (v .^ 2) ./ S
Z3 = Z[3]       # The load is at the third bus

# Define state vector

# Define parameter vector


# -----------------------------------------------------------------------------------------
# Section 3: Build ODEProblem
# -----------------------------------------------------------------------------------------

# Define function to be used by integration method 
function three_bus_network(du, u, p, t)
    ##### TODO: Write the system of equations  
    i12d, i12q, i13d, i13q, i23d, i23q, v1d, v1q, v2d, v2q, v3d, v3q = u
    R12, X12, B12, R13, X13, B13, R23, X23, B23, i1d, i1q, i2d, i2q, i3d, i3q = p

    # Algebraic Equations:

    # Line 1-2 system of equations:
    du[1] = (v1d - R12 * i12d + X12 * i12q) / X12
    du[2] = (v1q - R12 * i12q - X12 * i12d) / X12
    du[3] = (-B12 / 2) * v1q
    du[4] = (-B12 / 2) * v1d
    du[5] = (-B12 / 2) * v2q
    du[6] = (-B12 / 2) * v2d

    # Line 2-3 system of equations:
    du[7] = (v2d - R23 * i23d + X23 * i23q) / X23   
    du[8] = (v2q - R23 * i23q - X23 * i23d) / X23
    du[9] = (-B23 / 2) * v2q
    du[10] = (-B23 / 2) * v2d
    du[11] = (-B23 / 2) * v3q
    du[12] = (-B23 / 2) * v3d

    # Line 3-1 system of equations
    du[13] = (v3d - R31 * i13d + X13 * i13q) / X13
    du[14] = (v3q - R31 * i13q - X13 * i13d) / X13
    du[15] = (-B13 / 2) * v3q
    du[16] = (-B13 / 2) * v3d
    du[17] = (-B13 / 2) * v1q
    du[18] = (-B13 / 2) * v1d
end

# Build the mass matrix
M_diagonal = [] ##### TODO: Fill this with the coefficients of the derivatives 
M = Diagonal(M_diagonal)  # creates diagonal square matrix from a vector

# Build function 
f = ODEFunction(three_bus_network, mass_matrix=M)

# Define length of simulation
tspan = (0.0, 0.2)

# Build initial condition vector
u0 = [] ##### TODO: Fill this with initial conditions you calculated in Section 2

# Build parameter vector
p = [] ##### TODO: Fill this with parameters you calculated in Section 2

# Check initial condition consistency
# NOTE: We want this pre-perturbation initial condition to be an equilibrium point
du0 = ones(Float64, length(u0)); # return object for three_bus_network()
three_bus_network(du0, u0, p, 0); # evaluate function once
if norm(du0, Inf) > 1e-10
    throw("Residual norm of equilibrium IC is probably too big: $(norm(du0,Inf))")
else
    println("Residual norm of equilibrium IC looks good: $(norm(du0,Inf))")
end

# Build problem
prob = ODEProblem(f, u0, tspan, p)


# -----------------------------------------------------------------------------------------
# Section 4: Define perturbation and run simulation
# -----------------------------------------------------------------------------------------

# Define the set of times to apply a perturbation 
perturb_times = [0.01]

# Define the condition for which to apply a perturbation 
function condition(u, t, integrator)
    t in perturb_times
end

# Define the perturbation 

function affect!(integrator)
    integrator.p[999] *= 1.5 ##### TODO: Change 999 to index of the load impedance parameter
end

# Create a Callback function that represents the pertubation 
cb = DiscreteCallback(condition, affect!)

# Run simulation
sol = solve(prob, Rodas5P(), callback=cb, tstops=perturb_times)


# -----------------------------------------------------------------------------------------
# Section 5: Plot
# -----------------------------------------------------------------------------------------

# Plot (examples)
plot(sol,
    title="Three Bus Network: All Variables",
    xlabel="Time [s]",
)
plot(sol,
    idxs=(0, [1, 2]),
    title="Three Bus Network: Subset of Variables",
    xlabel="Time [s]",
    label=["label for u[1]" "label for u[2]"],
)
