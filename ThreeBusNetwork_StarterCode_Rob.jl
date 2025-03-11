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

# Find complex power
S = complex.(P, Q)

# Find complex voltage in ABC reference frame
V_abc = v .* ℯ .^ (im .* θ)

##### MOVE TO NETWORK DQ REFERENCE FRAME #####
# Find DQ voltage (Milano Eigenvalue Problems Eqn 1.44)
v_d = -v .* sin.(θ)             # Direct-axis component of bus voltage
v_q = v .* cos.(θ)              # Quadrature-axis component of bus voltage
V_dq = complex.(v_d, v_q)       # Complex bus voltage in network DQ reference frame

# Sanity check
V_test = V_dq .* ℯ .^ (-im * π / 2)
if norm(V_test - V_abc, Inf) > 1e-10
    throw("DQ voltage calculation is probably wrong. Difference between calculated and expected: $(norm(V_test - V_abc,Inf))")
else
    println("DQ voltage calculation looks good. Difference between calculated and expected: $(norm(V_test - V_abc,Inf))")
end

# Find complex current
I_dq = conj(S ./ V_dq)          # Complex bus current injections in network DQ reference frame
i_d = real(I_dq)                # Direct-axis component of bus current injections
i_q = imag(I_dq)                # Quadrature-axis component of bus current injections

# Sanity check
P_test = (v_d .* i_d) .+ (v_q .* i_q)
Q_test = (v_q .* i_d) .- (v_d .* i_q)
S_test = complex.(P_test, Q_test)
if norm(S_test - S, Inf) > 1e-10
    throw("DQ current calculation is probably wrong. Difference between calculated and expected apparent power: $(norm(S_test - S,Inf))")
else
    println("DQ current calculation looks good. Difference between calculated and expected apparent power: $(norm(S_test - S,Inf))")
end

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
i_1_d = i_d[1]
i_1_q = i_q[1]
i_2_d = i_d[2]
i_2_q = i_q[2]
i_3_d = i_d[3]
i_3_q = i_q[3]

# Define state vector elements (Assume steady-state)
v_1_d = v_d[1]                                      # Bus 1 d-axis terminal voltage             
v_1_q = v_q[1]                                      # Bus 1 q-axis terminal voltage             
v_2_d = v_d[2]                                      # Bus 2 d-axis terminal voltage             
v_2_q = v_q[2]                                      # Bus 2 q-axis terminal voltage             
v_3_d = v_d[3]                                      # Bus 3 d-axis terminal voltage             
v_3_q = v_q[3]                                      # Bus 3 q-axis terminal voltage  
i_b12_d_1 = -1 * (B_12 / 2) * v_1_q                 # Line 1-2 d-axis shunt current at bus 1
i_b12_q_1 = (B_12 / 2) * v_1_d                      # Line 1-2 q-axis shunt current at bus 1
i_b12_d_2 = -1 * (B_12 / 2) * v_2_q                 # Line 1-2 d-axis shunt current at bus 2
i_b12_q_2 = (B_12 / 2) * v_2_d                      # Line 1-2 q-axis shunt current at bus 2
i_b13_d_1 = -1 * (B_13 / 2) * v_1_q                 # Line 1-3 d-axis shunt current at bus 1
i_b13_q_1 = (B_13 / 2) * v_1_d                      # Line 1-3 q-axis shunt current at bus 1
i_b13_d_3 = -1 * (B_13 / 2) * v_3_q                 # Line 1-3 d-axis shunt current at bus 3
i_b13_q_3 = (B_13 / 2) * v_3_d                      # Line 1-3 q-axis shunt current at bus 3
i_b23_d_2 = -1 * (B_23 / 2) * v_2_q                 # Line 2-3 d-axis shunt current at bus 2
i_b23_q_2 = (B_23 / 2) * v_1_d                      # Line 2-3 q-axis shunt current at bus 2
i_b23_d_3 = -1 * (B_23 / 2) * v_3_q                 # Line 2-3 d-axis shunt current at bus 3
i_b23_q_3 = (B_23 / 2) * v_3_d                      # Line 2-3 q-axis shunt current at bus 3

# Still working on this – using the pseudoinverse for now because the A matrix is singular
A = [1 1 0; 1 0 -1; 0 1 1]
b = [i_1_d - (i_b12_d_1 + i_b13_d_1); (i_b12_d_2 + i_b23_d_2) - i_2_d; (i_b13_d_3 + i_b23_d_3) - i_3_d]
(i_12_d, i_13_d, i_23_d) = pinv(A) * b

A = [1 1 0; 1 0 -1; 0 1 1]
b = [i_1_q - (i_b12_q_1 + i_b13_q_1); (i_b12_q_2 + i_b23_q_2) - i_2_q; (i_b13_q_3 + i_b23_q_3) - i_3_q]
(i_12_q, i_13_q, i_23_q) = pinv(A) * b

# Sanity check
res_i = [
    i_1_d - i_12_d - i_13_d - i_b12_d_1 - i_b13_d_1;
    i_1_q - i_12_q - i_13_q - i_b12_q_1 - i_b13_q_1;
    i_2_d + i_12_d - i_23_d - i_b12_d_2 - i_b23_d_2;
    i_2_q + i_12_q - i_23_q - i_b12_q_2 - i_b23_q_2;
    i_3_d + i_23_d + i_13_d - i_b13_d_3 - i_b23_d_3;
    i_3_q + i_23_q + i_13_q - i_b13_q_3 - i_b23_q_3
]

if norm(res_i, Inf) > 1e-10
    @warn "DQ current calculation is probably wrong. Unaccounted for current: $(norm(res_i,Inf))"
else
    println("DQ current calculation looks good. Unaccounted for current: $(norm(res_i - S,Inf))")
end

# -----------------------------------------------------------------------------------------
# Section 3: Build ODEProblem
# -----------------------------------------------------------------------------------------

# Define function to be used by integration method 
function three_bus_network(du, u, p, t)
    i_12_d, i_12_q, i_13_d, i_13_q, i_23_d, i_23_q, v_1_d, v_1_q, v_2_d, v_2_q, v_3_d, v_3_q, i_b12_d_1, i_b12_q_1, i_b12_d_2, i_b12_q_2, i_b13_d_1, i_b13_q_1, i_b13_d_3, i_b13_q_3, i_b23_d_2, i_b23_q_2, i_b23_d_3, i_b23_q_3 = u
    R_12, X_12, B_12, R_13, X_13, B_13, R_23, X_23, B_23, i_1_d, i_1_q, i_2_d, i_2_q, i_3_d, i_3_q = p

    # Helper variables
    B1_eq = (B_12 + B_13)
    B2_eq = (B_12 + B_23)
    B3_eq = (B_13 + B_23)

    # Define equations in the same order as the u vector
    # Line currents (differential)
    du[1] = (v_1_d - v_2_d - R_12 * i_12_d + X_12 * i_12_q) / X_12          # d/dt (i_12_d)
    du[2] = (v_1_q - v_2_q - R_12 * i_12_q - X_12 * i_12_d) / X_12          # d/dt (i_12_q)
    du[3] = (v_1_d - v_3_d - R_13 * i_13_d + X_13 * i_13_q) / X_13          # d/dt (i_13_d)
    du[4] = (v_1_q - v_3_q - R_13 * i_13_q - X_13 * i_13_d) / X_13          # d/dt (i_13_q)
    du[5] = (v_2_d - v_3_d - R_23 * i_23_d + X_23 * i_23_q) / X_23          # d/dt (i_23_d)
    du[6] = (v_2_q - v_3_q - R_23 * i_23_q - X_23 * i_23_d) / X_23          # d/dt (i_23_q)

    # Bus voltages (differential)
    du[7] = (B1_eq / 2) * v_1_q + i_b12_d_1 + i_b13_d_1                     # d/dt (v_1_d)
    du[8] = -1 * (B1_eq / 2) * v_1_d + i_b12_q_1 + i_b13_q_1                # d/dt (v_1_q)
    du[9] = (B2_eq / 2) * v_2_q + i_b12_d_2 + i_b23_d_2                     # d/dt (v_2_d)
    du[10] = -1 * (B2_eq / 2) * v_2_d + i_b12_q_2 + i_b23_q_2               # d/dt (v_2_q)
    du[11] = (B3_eq / 2) * v_3_q + i_b23_d_3 + i_b13_d_3                    # d/dt (v_3_d)
    du[12] = -1 * (B3_eq / 2) * v_3_d + i_b23_q_3 + i_b13_q_3               # d/dt (v_3_q)

    # Shunt currents (algebraic)
    # Line 1-2
    du[13] = i_1_d - i_12_d - i_13_d - i_b12_d_1 - i_b13_d_1                 # d/dt (i_b12_d_1) = 0
    du[14] = i_1_q - i_12_q - i_13_q - i_b12_q_1 - i_b13_q_1                 # d/dt (i_b12_q_1) = 0
    du[15] = i_2_d + i_12_d - i_23_d - i_b12_d_2 - i_b23_d_2                 # d/dt (i_b12_d_2) = 0
    du[16] = i_2_q + i_12_q - i_23_q - i_b12_q_2 - i_b23_q_2                 # d/dt (i_b12_q_2) = 0
    # Line 1-3
    du[17] = i_1_d - i_12_d - i_13_d - i_b12_d_1 - i_b13_d_1                 # d/dt (i_b13_d_1) = 0
    du[18] = i_1_q - i_12_q - i_13_q - i_b12_q_1 - i_b13_q_1                 # d/dt (i_b13_q_1) = 0
    du[19] = i_3_d + i_23_d + i_13_d - i_b23_d_3 - i_b13_d_3                 # d/dt (i_b13_d_3) = 0
    du[20] = i_3_q + i_23_q + i_13_q - i_b23_q_3 - i_b13_q_3                 # d/dt (i_b13_q_3) = 0
    # Line 2-3
    du[21] = i_2_d + i_12_d - i_23_d - i_b12_d_2 - i_b23_d_2                 # d/dt (i_b23_d_2) = 0
    du[22] = i_2_q + i_12_q - i_23_q - i_b12_q_2 - i_b23_q_2                 # d/dt (i_b23_q_2) = 0
    du[23] = i_3_d + i_23_d + i_13_d - i_b23_d_3 - i_b13_d_3                 # d/dt (i_b23_d_3) = 0
    du[24] = i_3_q + i_23_q + i_13_q - i_b23_q_3 - i_b13_q_3                 # d/dt (i_b23_q_3) = 0
end

# Build the mass matrix
M_diagonal = [ones(12); zeros(12)]
M = Diagonal(M_diagonal)  # creates diagonal square matrix from a vector

# Build function 
f = ODEFunction(three_bus_network, mass_matrix=M)

# Define length of simulation
tspan = (0.0, 0.2)

# Build initial condition vector
# Order will be:
#u0= [i_12_d, i_12_q, i_13_d, i_13_q, i_23_d, i_23_q, v_1_d, v_1_q, v_2_d, v_2_q, v_3_d, v_3_q, i_b12_d_1, i_b12_q_1, i_b12_d_2, i_b12_q_2, i_b13_d_1, i_b13_q_1, i_b13_d_3, i_b13_q_3, i_b23_d_2, i_b23_q_2, i_b23_d_3, i_b23_q_3]
u0 = [i_12_d, i_12_q, i_13_d, i_13_q, i_23_d, i_23_q, v_1_d, v_1_q, v_2_d, v_2_q, v_3_d, v_3_q, i_b12_d_1, i_b12_q_1, i_b12_d_2, i_b12_q_2, i_b13_d_1, i_b13_q_1, i_b13_d_3, i_b13_q_3, i_b23_d_2, i_b23_q_2, i_b23_d_3, i_b23_q_3]


# Build parameter vector
# Order will be:
#p= [R_12, X_12, B_12, R_13, X_13, B_13, R_23, X_23, B_23, i_1_d, i_1_q, i_2_d, i_2_q, i_3_d, i_3_q]
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
