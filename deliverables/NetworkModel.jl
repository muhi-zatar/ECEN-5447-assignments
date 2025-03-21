# Defining the module
module NetworkModel

# Exporing the functions
export ThreeBusNetwork, initialize_network, update_network_states!, BUS_MACHINE_MODEL, sanity_check

using LinearAlgebra

# Definining some constants
const BUS_INFINITE_BUS = 1
const BUS_MACHINE_MODEL = 2
const BUS_LOAD = 3
const NUM_STATES = 16
const I_12_D_IDX = 1
const I_12_Q_IDX = 2
const I_13_D_IDX = 3
const I_13_Q_IDX = 4
const I_23_D_IDX = 5
const I_23_Q_IDX = 6
const V_1_D_IDX = 7
const V_1_Q_IDX = 8
const V_2_D_IDX = 9
const V_2_Q_IDX = 10
const V_3_D_IDX = 11
const V_3_Q_IDX = 12
const I_1_D_IDX = 13
const I_1_Q_IDX = 14
const I_3_D_IDX = 15
const I_3_Q_IDX = 16

# Helper functions
# not all of them are used.

function sanity_check(test_value, true_value, calculation_under_test::String, verbose=true)
    difference = norm(test_value .- true_value, Inf)
    if (difference > 1e-6)
        throw("$calculation_under_test calculation is probably wrong. Difference between calculated and expected: $difference")
    else
        if verbose
            println("$calculation_under_test calculation looks good. Difference between calculated and expected: $difference")
        end
    end
end

# RI to DQ transformation matrix
function ri_dq(delta)
    return [cos(delta) sin(delta); -sin(delta) cos(delta)]
end

# DQ to RI transformation matrix
function dq_ri(delta)
    return [cos(delta) -sin(delta); sin(delta) cos(delta)]
end

# Network structure for parameters
mutable struct ThreeBusNetwork
    # Network parameters
    R_12::Float64                                 # Line 1-2 resistance
    X_12::Float64                                 # Line 1-2 reactance
    B_1::Float64                                  # Lumped susceptance at Bus 1
    R_13::Float64                                 # Line 1-3 resistance
    X_13::Float64                                 # Line 1-3 reactance
    B_3::Float64                                  # Lumped susceptance at Bus 3
    R_23::Float64                                 # Line 2-3 resistance
    X_23::Float64                                 # Line 2-3 reactance
    B_2::Float64                                  # Lumped susceptance at Bus 2
    X_IB::Float64                                 # Reactance between network and infinite bus EMF
    E_IB_D::Float64                               # d-axis EMF behind reactance of the infinite bus
    E_IB_Q::Float64                               # q-axis EMF behind reactance of the infinite bus
    Z_L::Complex{Float64}                         # Impedance of the load at Bus 3
    #M::AbstractArray{Float64}                     # The coefficients that go on the diagonal of the mass matrix

    # Constructor with default values
    function ThreeBusNetwork(;
        R_12=0.01,
        X_12=0.12,
        B_1=0.05,
        R_13=0.01,
        X_13=0.12,
        B_3=0.05,
        R_23=0.01,
        X_23=0.12,
        B_2=0.05,
        X_IB=0.1,
        E_IB_D=0.0,                               # To be filled in from power flow solution
        E_IB_Q=0.0,                               # To be filled in from power flow solution
        Z_L=0.0 + im * 0.0,                       # To be filled in from power flow solution
    )
        return new(R_12, X_12, B_1, R_13, X_13, B_3, R_23, X_23, B_2, X_IB, E_IB_D, E_IB_Q, Z_L)
    end # constructor function
end # mutable struct

```
The module has two functions:
Initializing states, and updating states
This will help us in scalability
```

# Initialize network states based on power flow solution
```
This function will initialize the states of the network at steady-state.
Two important notes when working with the network:
    1 -- We work in the DQ reference frame of the network, which is different than the DQ reference
       frame of the machine. The network DQ reference frame rotates at a constant angular velocity,
       which is assumed to be a constant 60Hz.
    2 -- We work on every bus in the network, which means that we want the power flow results to be
       vectors of voltages, angles, etc., rather than single values as used by the generator 
       components.
```
function initialize_network(network::ThreeBusNetwork, V_m::Vector{Float64}, θ::Vector{Float64}, P::Vector{Float64}, Q::Vector{Float64})
    # Prepare state vector for return
    states = zeros(Float64, NUM_STATES)

    # Reconstruct complex power and terminal voltage
    S = complex.(P, Q)
    V_terminal = V_m .* ℯ .^ (im .* θ)


    print("\n\nInitial positive sequence voltage at the bus: $(V_m[BUS_MACHINE_MODEL]) ∠ $(rad2deg(θ[BUS_MACHINE_MODEL]))\n\n")

    # Construct initial current to pass forward
    I_RI = conj.(S ./ V_terminal)
    I_2_RI = I_RI[BUS_MACHINE_MODEL]
    I_2_m = abs(I_2_RI)
    I_2_ang = angle(I_2_RI)
    print("\n\nInitial positive sequence current at the bus: $(I_2_m) ∠ $(rad2deg(I_2_ang))\n\n")

    ##### MOVE TO NETWORK DQ REFERENCE FRAME #####
    # Find DQ voltage (Milano Eigenvalue Problems Eqn 1.44)
    v_d = -V_m .* sin.(θ)             # Direct-axis component of bus voltage
    v_q = V_m .* cos.(θ)              # Quadrature-axis component of bus voltage
    V_dq = complex.(v_d, v_q)         # Complex bus voltages in network DQ reference frame

    # Sanity check
    V_test = V_dq .* ℯ .^ (-im * π / 2)
    sanity_check(V_test, V_terminal, "DQ voltage")

    # Find complex current
    I_dq = conj(S ./ V_dq)            # Complex bus current injections in network DQ reference frame
    i_d = real(I_dq)                  # Direct-axis component of bus current injections
    i_q = imag(I_dq)                  # Quadrature-axis component of bus current injections
    I_test = I_dq .* ℯ .^ (-im * π / 2)
    sanity_check(I_test, I_RI, "DQ Current")

    # Bus 2 (our model) injection current in network reference frame
    i_2_d = i_d[BUS_MACHINE_MODEL]
    i_2_q = i_q[BUS_MACHINE_MODEL]

    # Sanity check
    P_test = (v_d .* i_d) .+ (v_q .* i_q)
    Q_test = (v_q .* i_d) .- (v_d .* i_q)
    S_test = complex.(P_test, Q_test)
    sanity_check(S_test, S, "Power")

    # Find load impedance
    Z_dq = (abs.(V_dq) .^ 2) ./ conj.(S)
    network.Z_L = Z_dq[BUS_LOAD]             # The load is at the third bus

    # Define state vector elements (Assume steady-state)
    states[I_1_D_IDX] = i_d[1]                                      # Bus 1 d-axis injected current
    states[I_1_Q_IDX] = i_q[1]                                      # Bus 1 q-axis injected current
    states[I_3_D_IDX] = i_d[3]                                      # Bus 3 d-axis injected current
    states[I_3_Q_IDX] = i_q[3]                                      # Bus 3 q-axis injected current
    states[V_1_D_IDX] = v_d[1]                                      # Bus 1 d-axis terminal voltage             
    states[V_1_Q_IDX] = v_q[1]                                      # Bus 1 q-axis terminal voltage             
    states[V_2_D_IDX] = v_d[2]                                      # Bus 2 d-axis terminal voltage             
    states[V_2_Q_IDX] = v_q[2]                                      # Bus 2 q-axis terminal voltage             
    states[V_3_D_IDX] = v_d[3]                                      # Bus 3 d-axis terminal voltage             
    states[V_3_Q_IDX] = v_q[3]                                      # Bus 3 q-axis terminal voltage  
    i_b1_d = -1 * network.B_1 * v_q[1]                  # Bus 1 d-axis shunt current
    i_b1_q = network.B_1 * v_d[1]                       # Bus 1 q-axis shunt current
    i_b2_d = -1 * network.B_2 * v_q[2]                  # Bus 2 d-axis shunt current
    i_b2_q = network.B_2 * v_d[2]                       # Bus 2 q-axis shunt current
    i_b3_d = -1 * network.B_3 * v_q[3]                  # Bus 3 d-axis shunt current
    i_b3_q = network.B_3 * v_d[3]                       # Bus 3 q-axis shunt current

    # Voltage behind reactance
    # TODO: If eliminating mass matrix doesn't work, look into this
    network.E_IB_D = v_d[BUS_INFINITE_BUS] - i_q[BUS_INFINITE_BUS] * network.X_IB
    network.E_IB_Q = v_q[BUS_INFINITE_BUS] + i_d[BUS_INFINITE_BUS] * network.X_IB

    # Sanity check on voltages behind reactance
    E1 = complex(network.E_IB_D, network.E_IB_Q)
    if abs(E1) < V_m[BUS_INFINITE_BUS]
        throw("Voltage behind reactance is less than terminal voltage at Bus 1! Check currents...")
    end

    # Line Currents
    A = [network.R_12 -1*network.X_12 0 0 0 0;
        network.X_12 network.R_12 0 0 0 0;
        0 0 network.R_13 -1*network.X_13 0 0;
        0 0 network.X_13 network.R_13 0 0;
        0 0 0 0 network.R_23 -1*network.X_23;
        0 0 0 0 network.X_23 network.R_23]
    b = [v_d[1] - v_d[2];
        v_q[1] - v_q[2];
        v_d[1] - v_d[3];
        v_q[1] - v_q[3];
        v_d[2] - v_d[3];
        v_q[2] - v_q[3]]
    line_currents = A \ b
    states[I_12_D_IDX] = line_currents[1]                           # Line 1-2 d-axis current
    states[I_12_Q_IDX] = line_currents[2]                           # Line 1-2 q-axis current
    states[I_13_D_IDX] = line_currents[3]                           # Line 1-3 d-axis current
    states[I_13_Q_IDX] = line_currents[4]                           # Line 1-3 q-axis current
    states[I_23_D_IDX] = line_currents[5]                           # Line 2-3 d-axis current
    states[I_23_Q_IDX] = line_currents[6]                           # Line 2-3 q-axis current

    # Sanity check on line currents (this is just a KVL on the calculated initial values)
    res_i = [
        states[I_1_D_IDX] - states[I_12_D_IDX] - states[I_13_D_IDX] - i_b1_d;
        states[I_1_Q_IDX] - states[I_12_Q_IDX] - states[I_13_Q_IDX] - i_b1_q;
        i_2_d + states[I_12_D_IDX] - states[I_23_D_IDX] - i_b2_d;
        i_2_q + states[I_12_Q_IDX] - states[I_23_Q_IDX] - i_b2_q;
        states[I_3_D_IDX] + states[I_23_D_IDX] + states[I_13_D_IDX] - i_b3_d;
        states[I_3_Q_IDX] + states[I_23_Q_IDX] + states[I_13_Q_IDX] - i_b3_q
    ]

    sanity_check(res_i, zeros(6), "Line current")

    # # Populate mass matrix
    # M_diagonal = zeros(Float64, NUM_STATES)
    # M_diagonal[I_12_D_IDX] = network.X_12
    # M_diagonal[I_12_Q_IDX] = network.X_12
    # M_diagonal[I_13_D_IDX] = network.X_13
    # M_diagonal[I_13_Q_IDX] = network.X_13
    # M_diagonal[I_23_D_IDX] = network.X_23
    # M_diagonal[I_23_Q_IDX] = network.X_23
    # M_diagonal[V_1_D_IDX] = network.B_1
    # M_diagonal[V_1_Q_IDX] = network.B_1
    # M_diagonal[V_2_D_IDX] = network.B_2
    # M_diagonal[V_2_Q_IDX] = network.B_2
    # M_diagonal[V_3_D_IDX] = network.B_3
    # M_diagonal[V_3_Q_IDX] = network.B_3
    # M_diagonal[I_1_D_IDX] = network.X_IB
    # M_diagonal[I_1_Q_IDX] = network.X_IB
    # M_diagonal[I_3_D_IDX] = imag(network.Z_L)
    # M_diagonal[I_3_Q_IDX] = imag(network.Z_L)
    # M_diagonal[I_B1_D_IDX:I_B3_Q_IDX] .= 0.0
    # network.M = M_diagonal

    # Return initial states
    return states, i_2_d, i_2_q
end

function update_network_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    S_bus::Complex{Float64},
    network::ThreeBusNetwork
)
    # Extract network states
    i_12_d = states[I_12_D_IDX]
    i_12_q = states[I_12_Q_IDX]
    i_13_d = states[I_13_D_IDX]
    i_13_q = states[I_13_Q_IDX]
    i_23_d = states[I_23_D_IDX]
    i_23_q = states[I_23_Q_IDX]
    v_1_d = states[V_1_D_IDX]
    v_1_q = states[V_1_Q_IDX]
    v_2_d = states[V_2_D_IDX]
    v_2_q = states[V_2_Q_IDX]
    v_3_d = states[V_3_D_IDX]
    v_3_q = states[V_3_Q_IDX]
    i_1_d = states[I_1_D_IDX]
    i_1_q = states[I_1_Q_IDX]
    i_3_d = states[I_3_D_IDX]
    i_3_q = states[I_3_Q_IDX]
    # i_b1_d = states[I_B1_D_IDX]
    # i_b1_q = states[I_B1_Q_IDX]
    # i_b2_d = states[I_B2_D_IDX]
    # i_b2_q = states[I_B2_Q_IDX]
    # i_b3_d = states[I_B3_D_IDX]
    # i_b3_q = states[I_B3_Q_IDX]

    ## Calculate network DQ current from bus power and voltage
    V_dq = complex(v_2_d, v_2_q)
    I_dq = conj(S_bus ./ V_dq)            # Complex bus current injections in network DQ reference frame
    i_2_d = real(I_dq)                    # Direct-axis component of bus current injections
    i_2_q = imag(I_dq)                    # Quadrature-axis component of bus current injections

    # Shunt currents (algebraic)
    # Bus 1
    i_b1_d = i_1_d - i_12_d - i_13_d                                                # d/dt (i_b1_d) = 0
    i_b1_q = i_1_q - i_12_q - i_13_q                                                # d/dt (i_b1_q) = 0
    # Bus 2
    i_b2_d = i_2_d + i_12_d - i_23_d                                                # d/dt (i_b2_d) = 0
    i_b2_q = i_2_q + i_12_q - i_23_q                                                # d/dt (i_b2_q) = 0
    # Bus 3
    i_b3_d = i_3_d + i_23_d + i_13_d                                                 # d/dt (i_b3_d) = 0
    i_b3_q = i_3_q + i_23_q + i_13_q                                                 # d/dt (i_b3_q) = 0

    # Unpack real and imaginary components of the load for convenience
    R_load = real(network.Z_L)
    X_load = imag(network.Z_L)

    # Define derivatives of all states (differential and algebraic)
    # Line currents (differential)
    derivatives[I_12_D_IDX] = (v_1_d - v_2_d - network.R_12 * i_12_d + network.X_12 * i_12_q) / network.X_12                 # d/dt (i_12_d) != 0
    derivatives[I_12_Q_IDX] = (v_1_q - v_2_q - network.R_12 * i_12_q - network.X_12 * i_12_d) / network.X_12                 # d/dt (i_12_q) != 0
    derivatives[I_13_D_IDX] = (v_1_d - v_3_d - network.R_13 * i_13_d + network.X_13 * i_13_q) / network.X_13                 # d/dt (i_13_d) != 0
    derivatives[I_13_Q_IDX] = (v_1_q - v_3_q - network.R_13 * i_13_q - network.X_13 * i_13_d) / network.X_13                 # d/dt (i_13_q) != 0
    derivatives[I_23_D_IDX] = (v_2_d - v_3_d - network.R_23 * i_23_d + network.X_23 * i_23_q) / network.X_23                 # d/dt (i_23_d) != 0
    derivatives[I_23_Q_IDX] = (v_2_q - v_3_q - network.R_23 * i_23_q - network.X_23 * i_23_d) / network.X_23                 # d/dt (i_23_q) != 0

    # Bus voltages (differential)
    derivatives[V_1_D_IDX] = (i_b1_d + network.B_1 * v_1_q) / network.B_1                                                     # d/dt (v_1_d) != 0
    derivatives[V_1_Q_IDX] = (i_b1_q - network.B_1 * v_1_d) / network.B_1                                                     # d/dt (v_1_q) != 0
    derivatives[V_2_D_IDX] = (i_b2_d + network.B_2 * v_2_q) / network.B_2                                                     # d/dt (v_2_d) != 0
    derivatives[V_2_Q_IDX] = (i_b2_q - network.B_2 * v_2_d) / network.B_2                                                     # d/dt (v_2_q) != 0
    derivatives[V_3_D_IDX] = (i_b3_d + network.B_3 * v_3_q) / network.B_3                                                     # d/dt (v_3_d) != 0
    derivatives[V_3_Q_IDX] = (i_b3_q - network.B_3 * v_3_d) / network.B_3                                                     # d/dt (v_3_q) != 0

    # Current injections
    # Bus 1 (Infinite Bus) (differential) –– From Milano's Eigenvalue Problems book, Eq. 1.48
    derivatives[I_1_D_IDX] = ((network.E_IB_D - v_1_d) + network.X_IB * i_1_q) / network.X_IB                                  # d/dt (i_1_d) != 0
    derivatives[I_1_Q_IDX] = ((network.E_IB_Q - v_1_q) - network.X_IB * i_1_d) / network.X_IB                                  # d/dt (i_1_q) != 0
    # Bus 3 (Load) (differential) –– From Milano's Eigenvalue Problems book, Eq. 1.48
    derivatives[I_3_D_IDX] = (v_3_d - R_load * i_3_d + X_load * i_3_q) / imag(network.Z_L)                                          # d/dt (i_3_d) != 0
    derivatives[I_3_Q_IDX] = (v_3_q - R_load * i_3_q - X_load * i_3_d) / imag(network.Z_L)                                          # d/dt (i_3_q) != 0

    # Compute apparent power to return (Milano Eigenvalue Problems Eq. 1.42)
    P_terminal = v_2_d * i_2_d + v_2_q * i_2_q
    Q_terminal = v_2_q * i_2_d - v_2_d * i_2_q
    S_terminal = complex(P_terminal, Q_terminal)

    # Compute positive sequence terminal voltage to return (Milano Eigenvalue Problems Eq 1.43)
    V_terminal = (v_2_d + im * v_2_q) * exp(-im * π / 2)

    # Compute positive sequence terminal injected current to return (can check this against Milano Eigenvalue Problems Eq 1.41)
    I_terminal = conj(S_terminal / V_terminal)
    I_terminal_test = (i_2_d + im * i_2_q) * exp(-im * π / 2)
    sanity_check(I_terminal_test, I_terminal, "Bus 2 Current Injection", false)

    # Return all values (used for debugging – we can return fewer values after we figure out what's wrong)
    return V_terminal, S_terminal, I_terminal, i_2_d, i_2_q
end
end # module