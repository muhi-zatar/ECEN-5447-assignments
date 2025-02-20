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

# Define some arbitrary generator reactances
x1 = 1.0
x2 = 1.0

# Find load impedance
Z = (v .^ 2) ./ S
Z3 = Z[3]       # The load is at the third bus

# Find Lines
line_1_2 = get_component(Line, sys, "BUS 1-BUS 2-i_1")
line_1_3 = get_component(Line, sys, "BUS 1-BUS 3-i_1")
line_2_3 = get_component(Line, sys, "BUS 2-BUS 3-i_1")

# Define state vector elements
i_12_d = 0
i_12_q = 0
i_13_d = 0
i_13_q = 0
i_23_d = 0
i_23_q = 0
v_1_d = real(V)[1]
v_1_q = imag(V)[1]
v_2_d = real(V)[2]
v_2_q = imag(V)[1]
v_3_d = real(V)[3]
v_3_q = imag(V)[1]
i_b12_d_1 = 0
i_b12_q_1 = 0
i_b12_d_2 = 0
i_b12_q_2 = 0
i_b13_d_1 = 0
i_b13_q_1 = 0
i_b13_d_3 = 0
i_b13_q_3 = 0
i_b23_d_2 = 0
i_b23_q_2 = 0
i_b23_d_3 = 0
i_b23_q_3 = 0

# Define parameter vector elements
R_12 = get_r(line_1_2)
X_12 = get_x(line_1_2)
B_12 = get_b(line_1_2)[1] * 2
R_13 = get_r(line_1_3)
X_13 = get_x(line_1_3)
B_13 = get_b(line_1_3)[1] * 2
R_23 = get_r(line_2_3)
X_23 = get_x(line_2_3)
B_23 = get_b(line_2_3)[1] * 2
i_1_d = real(I)[1]
i_1_q = imag(I)[1]
i_2_d = real(I)[2]
i_2_q = imag(I)[2]
i_3_d = real(I)[3]
i_3_q = imag(I)[3]

# -----------------------------------------------------------------------------------------
# Section 3: Build ODEProblem
# -----------------------------------------------------------------------------------------

# Define function to be used by integration method 
function three_bus_network(du, u, p, t)
    ##### TODO: Write the system of equations 
end

# Build the mass matrix
M_diagonal = [] ##### TODO: Fill this with the coefficients of the derivatives 
M = Diagonal(M_diagonal)  # creates diagonal square matrix from a vector

# Build function 
f = ODEFunction(three_bus_network, mass_matrix=M)

# Define length of simulation
tspan = (0.0, 0.2)

# Build initial condition vector
# Order will be:
# [i_12_d, i_12_q, i_13_d, i_13_q, i_23_d, i_23_q, v_1_d, v_1_q, v_2_d, v_2_q, v_3_d, v_3_q, i_b12_d_1, i_b12_q_1, i_b12_d_2, i_b12_q_2, i_b13_d_1, i_b13_q_1, i_b13_d_3, i_b13_q_3, i_b23_d_2, i_b23_q_2, i_b23_d_3, i_b23_q_3]
##### TODO: Fill this with initial conditions you calculated in Section 2
u0 = [i_12_d, i_12_q, i_13_d, i_13_q, i_23_d, i_23_q, v_1_d, v_1_q, v_2_d, v_2_q, v_3_d, v_3_q, i_b12_d_1, i_b12_q_1, i_b12_d_2, i_b12_q_2, i_b13_d_1, i_b13_q_1, i_b13_d_3, i_b13_q_3, i_b23_d_2, i_b23_q_2, i_b23_d_3, i_b23_q_3]


# Build parameter vector
# Order will be:
# [R_12, X_12, B_12, R_13, X_13, B_13, R_23, X_23, B_23, i_1_d, i_1_q, i_2_d, i_2_q, i_3_d, i_3_q]
##### TODO: Fill this with parameters you calculated in Section 2
p = [R_12, X_12, B_12, R_13, X_13, B_13, R_23, X_23, B_23, i_1_d, i_1_q, i_2_d, i_2_q, i_3_d, i_3_q]

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
