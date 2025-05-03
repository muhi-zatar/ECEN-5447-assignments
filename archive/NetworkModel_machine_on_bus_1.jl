# Defining the module
module NetworkModel_machine_on_bus_1

# Exporing the functions
export ThreeBusNetwork, initialize_network, update_network_states!, BUS_MACHINE_MODEL, sanity_check

using LinearAlgebra

# Definining some constants
const BUS_INFINITE_BUS = 2
const BUS_MACHINE_MODEL = 1
const BUS_LOAD = 3
const NUM_STATES = 22
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
const I_2_D_IDX = 13
const I_2_Q_IDX = 14
const I_3_D_IDX = 15
const I_3_Q_IDX = 16
const I_B1_D_IDX = 17
const I_B1_Q_IDX = 18
const I_B2_D_IDX = 19
const I_B2_Q_IDX = 20
const I_B3_D_IDX = 21
const I_B3_Q_IDX = 22

# Define common constants and variables
# Matrix transformation functions
# These are used to convert between rectangular and polar coordinates
# Same as power system package to avoid confusion in this part.


# Helper functions
# not all of them are used.

function sanity_check(test_value, true_value, calculation_under_test::String, verbose=true)
    distance = norm(test_value .- true_value, Inf)
    if (distance > 1e-2)
        difference = test_value .- true_value
        println(difference)
        throw("$calculation_under_test calculation is probably wrong. Residual is: $distance")
    else
        if verbose
            println("$calculation_under_test calculation looks good. Residual is: $distance")
        end
    end
end

function ri_dq(δ::T) where {T<:Number}
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return Float64[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

function dq_ri(δ::T) where {T<:Number}
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return Float64[
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

function dq_transformer(voltage, current, voltage_angle, current_angle)
    # First, calculate with ri_dq function
    v_dq = ri_dq(voltage_angle) * voltage
    i_dq = ri_dq(current_angle) * current

    return v_dq, i_dq
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
    M::AbstractArray{Float64}                     # The coefficients that go on the diagonal of the mass matrix

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
        M=zeros(Float64, NUM_STATES, NUM_STATES)  # To be filled in during initialization
    )
        return new(R_12, X_12, B_1, R_13, X_13, B_3, R_23, X_23, B_2, X_IB, E_IB_D, E_IB_Q, Z_L, M)
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

    # Reconstruct complex voltage and current at the terminal to calculate DQ
    S = complex.(P, Q)
    V_terminal = V_m .* ℯ .^ (im .* θ)
    print("\n\nInitial positive sequence voltage at the bus: $(V_terminal[2])\n\n")
    I_terminal = conj.(S ./ V_terminal)
    print("\n\nInitial positive sequence current at the bus: $(I_terminal[2])\n\n")

    ##### MOVE TO NETWORK DQ REFERENCE FRAME #####
    # Find DQ voltage (Milano Eigenvalue Problems Eqn 1.44)
    I_dq = zeros(Complex{Float64}, 3)
    V_dq = zeros(Complex{Float64}, 3)
    for i in 1:3
        V_RI = [real(V_terminal[i]); imag(V_terminal[i])]
        I_RI = [real(I_terminal[i]); imag(I_terminal[i])]
        V_DQ, I_DQ = dq_transformer(V_RI, I_RI, 0.0, 0.0)
        V_dq[i] = complex(V_DQ[1], V_DQ[2])
        I_dq[i] = complex(I_DQ[1], I_DQ[2])
    end
    i_d = real.(I_dq)
    i_q = imag.(I_dq)
    v_d = real.(V_dq)
    v_q = imag.(V_dq)

    # Bus 1 (our model) injection current in network reference frame
    i_1_d = i_d[BUS_MACHINE_MODEL]
    i_1_q = i_q[BUS_MACHINE_MODEL]

    # Sanity check
    P_test = (v_d .* i_d) .+ (v_q .* i_q)
    Q_test = (v_q .* i_d) .- (v_d .* i_q)
    S_test = complex.(P_test, Q_test)
    sanity_check(S_test, S, "Power")

    # Find load impedance
    Z_dq = (abs.(V_dq) .^ 2) ./ conj.(S)
    network.Z_L = Z_dq[BUS_LOAD]             # The load is at the third bus

    # Define state vector elements (Assume steady-state)
    states[I_2_D_IDX] = i_d[2]                                      # Bus 2 d-axis injected current
    states[I_2_Q_IDX] = i_q[2]                                      # Bus 2 q-axis injected current
    states[I_3_D_IDX] = i_d[3]                                      # Bus 3 d-axis injected current
    states[I_3_Q_IDX] = i_q[3]                                      # Bus 3 q-axis injected current
    states[V_1_D_IDX] = v_d[1]                                      # Bus 1 d-axis terminal voltage             
    states[V_1_Q_IDX] = v_q[1]                                      # Bus 1 q-axis terminal voltage             
    states[V_2_D_IDX] = v_d[2]                                      # Bus 2 d-axis terminal voltage             
    states[V_2_Q_IDX] = v_q[2]                                      # Bus 2 q-axis terminal voltage             
    states[V_3_D_IDX] = v_d[3]                                      # Bus 3 d-axis terminal voltage             
    states[V_3_Q_IDX] = v_q[3]                                      # Bus 3 q-axis terminal voltage  
    states[I_B1_D_IDX] = -1 * network.B_1 * v_q[1]                  # Bus 1 d-axis shunt current
    states[I_B1_Q_IDX] = network.B_1 * v_d[1]                       # Bus 1 q-axis shunt current
    states[I_B2_D_IDX] = -1 * network.B_2 * v_q[2]                  # Bus 2 d-axis shunt current
    states[I_B2_Q_IDX] = network.B_2 * v_d[2]                       # Bus 2 q-axis shunt current
    states[I_B3_D_IDX] = -1 * network.B_3 * v_q[3]                  # Bus 3 d-axis shunt current
    states[I_B3_Q_IDX] = network.B_3 * v_d[3]                       # Bus 3 q-axis shunt current

    # Voltage behind reactance
    network.E_IB_D = v_d[BUS_INFINITE_BUS] - i_q[BUS_INFINITE_BUS] * network.X_IB
    network.E_IB_Q = v_q[BUS_INFINITE_BUS] + i_d[BUS_INFINITE_BUS] * network.X_IB

    # Sanity check on voltages behind reactance
    E2 = complex(network.E_IB_D, network.E_IB_Q)
    if abs(E2) < V_m[BUS_INFINITE_BUS]
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

    # DEBUGGING – PRINT LINE CURRENTS
    I_12 = dq_ri(0) * [line_currents[1]; line_currents[2]]
    I_13 = dq_ri(0) * [line_currents[3]; line_currents[4]]
    I_23 = dq_ri(0) * [line_currents[5]; line_currents[6]]
    println("I_12 = $I_12")
    println("I_13 = $I_13")
    println("I_23 = $I_23")

    # Sanity check on line currents (this is just a KVL on the calculated initial values)
    res_i = [
        i_1_d - states[I_12_D_IDX] - states[I_13_D_IDX] - states[I_B1_D_IDX];
        i_1_q - states[I_12_Q_IDX] - states[I_13_Q_IDX] - states[I_B1_Q_IDX];
        states[I_2_D_IDX] + states[I_12_D_IDX] - states[I_23_D_IDX] - states[I_B2_D_IDX];
        states[I_2_Q_IDX] + states[I_12_Q_IDX] - states[I_23_Q_IDX] - states[I_B2_Q_IDX];
        states[I_3_D_IDX] + states[I_23_D_IDX] + states[I_13_D_IDX] - states[I_B3_D_IDX];
        states[I_3_Q_IDX] + states[I_23_Q_IDX] + states[I_13_Q_IDX] - states[I_B3_Q_IDX]
    ]

    sanity_check(res_i, zeros(6), "Line current")

    # Populate mass matrix
    ω_0 = 2.0 * π * 60.0
    M_diagonal = zeros(Float64, NUM_STATES)
    M_diagonal[I_12_D_IDX] = network.X_12 / ω_0
    M_diagonal[I_12_Q_IDX] = network.X_12 / ω_0
    M_diagonal[I_13_D_IDX] = network.X_13 / ω_0
    M_diagonal[I_13_Q_IDX] = network.X_13 / ω_0
    M_diagonal[I_23_D_IDX] = network.X_23 / ω_0
    M_diagonal[I_23_Q_IDX] = network.X_23 / ω_0
    M_diagonal[V_1_D_IDX] = network.B_1 / ω_0
    M_diagonal[V_1_Q_IDX] = network.B_1 / ω_0
    M_diagonal[V_2_D_IDX] = network.B_2 / ω_0
    M_diagonal[V_2_Q_IDX] = network.B_2 / ω_0
    M_diagonal[V_3_D_IDX] = network.B_3 / ω_0
    M_diagonal[V_3_Q_IDX] = network.B_3 / ω_0
    M_diagonal[I_2_D_IDX] = 0.0
    M_diagonal[I_2_Q_IDX] = 0.0
    M_diagonal[I_3_D_IDX] = 0.0
    M_diagonal[I_3_Q_IDX] = 0.0
    M_diagonal[I_B1_D_IDX:I_B3_Q_IDX] .= 0.0
    network.M = M_diagonal

    # Return initial states
    return states, i_1_d, i_1_q
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
    i_2_d = states[I_2_D_IDX]
    i_2_q = states[I_2_Q_IDX]
    i_3_d = states[I_3_D_IDX]
    i_3_q = states[I_3_Q_IDX]
    i_b1_d = states[I_B1_D_IDX]
    i_b1_q = states[I_B1_Q_IDX]
    i_b2_d = states[I_B2_D_IDX]
    i_b2_q = states[I_B2_Q_IDX]
    i_b3_d = states[I_B3_D_IDX]
    i_b3_q = states[I_B3_Q_IDX]

    ## Calculate network DQ current from bus power and voltage
    V_dq = complex(v_1_d, v_1_q)
    I_dq = conj(S_bus ./ V_dq)            # Complex bus current injections in network DQ reference frame
    # println("v_1_d: $v_1_d")
    # println("v_1_q: $v_1_q")
    # println("S (Bus 1): $S_bus")
    i_1_d = real(I_dq)                    # Direct-axis component of bus current injections
    i_1_q = imag(I_dq)                    # Quadrature-axis component of bus current injections

    # Debugging
    # println("v_2_d: $v_2_d")
    # println("v_2_q: $v_2_q")
    P_terminal_2 = v_2_d * i_2_d + v_2_q * i_2_q
    Q_terminal_2 = v_2_q * i_2_d - v_2_d * i_2_q
    S_terminal_2 = complex(P_terminal_2, Q_terminal_2)
    # println("S (Bus 2): $S_terminal_2")

    # Unpack real and imaginary components of the load for convenience
    R_load = real(network.Z_L)
    X_load = imag(network.Z_L)

    # Define derivatives of all states (differential and algebraic)
    # Line currents (differential)
    derivatives[I_12_D_IDX] = (v_1_d - v_2_d - network.R_12 * i_12_d + network.X_12 * i_12_q)                 # d/dt (i_12_d) != 0
    derivatives[I_12_Q_IDX] = (v_1_q - v_2_q - network.R_12 * i_12_q - network.X_12 * i_12_d)                 # d/dt (i_12_q) != 0
    derivatives[I_13_D_IDX] = (v_1_d - v_3_d - network.R_13 * i_13_d + network.X_13 * i_13_q)                 # d/dt (i_13_d) != 0
    derivatives[I_13_Q_IDX] = (v_1_q - v_3_q - network.R_13 * i_13_q - network.X_13 * i_13_d)                 # d/dt (i_13_q) != 0
    derivatives[I_23_D_IDX] = (v_2_d - v_3_d - network.R_23 * i_23_d + network.X_23 * i_23_q)                 # d/dt (i_23_d) != 0
    derivatives[I_23_Q_IDX] = (v_2_q - v_3_q - network.R_23 * i_23_q - network.X_23 * i_23_d)                 # d/dt (i_23_q) != 0

    # Bus voltages (differential)
    derivatives[V_1_D_IDX] = (i_b1_d + network.B_1 * v_1_q)                                                     # d/dt (v_1_d) != 0
    derivatives[V_1_Q_IDX] = (i_b1_q - network.B_1 * v_1_d)                                                     # d/dt (v_1_q) != 0
    derivatives[V_2_D_IDX] = (i_b2_d + network.B_2 * v_2_q)                                                     # d/dt (v_2_d) != 0
    derivatives[V_2_Q_IDX] = (i_b2_q - network.B_2 * v_2_d)                                                     # d/dt (v_2_q) != 0
    derivatives[V_3_D_IDX] = (i_b3_d + network.B_3 * v_3_q)                                                     # d/dt (v_3_d) != 0
    derivatives[V_3_Q_IDX] = (i_b3_q - network.B_3 * v_3_d)                                                     # d/dt (v_3_q) != 0

    # Current injections
    # Bus 1 (Infinite Bus) (algebraic) –– From Milano's Eigenvalue Problems book, Eq. 1.48
    derivatives[I_2_D_IDX] = ((network.E_IB_D - v_2_d) + network.X_IB * i_2_q)                                  # d/dt (i_2_d) = 0
    derivatives[I_2_Q_IDX] = ((network.E_IB_Q - v_2_q) - network.X_IB * i_2_d)                                  # d/dt (i_2_q) = 0
    # Bus 3 (Load) (algebraic) –– From Milano's Eigenvalue Problems book, Eq. 1.48
    derivatives[I_3_D_IDX] = (v_3_d - R_load * i_3_d + X_load * i_3_q)                                          # d/dt (i_3_d) = 0
    derivatives[I_3_Q_IDX] = (v_3_q - R_load * i_3_q - X_load * i_3_d)                                          # d/dt (i_3_q) = 0

    # Shunt currents (algebraic)
    # Bus 1
    derivatives[I_B1_D_IDX] = i_1_d - i_12_d - i_13_d - i_b1_d                                                # d/dt (i_b1_d) = 0
    derivatives[I_B1_Q_IDX] = i_1_q - i_12_q - i_13_q - i_b1_q                                                # d/dt (i_b1_q) = 0
    # println("I_B1_D = $i_b1_d")
    # println("I_B1_Q = $i_b1_q")
    # println("I_12_D = $i_12_d")
    # println("I_12_Q = $i_12_q")
    # println("I_13_D = $i_13_d")
    # println("I_13_Q = $i_13_q")
    # println("I_1_D = $i_1_d")
    # println("I_1_Q = $i_1_q")
    # exit()
    # Bus 2
    derivatives[I_B2_D_IDX] = i_2_d + i_12_d - i_23_d - i_b2_d                                                # d/dt (i_b2_d) = 0
    derivatives[I_B2_Q_IDX] = i_2_q + i_12_q - i_23_q - i_b2_q                                                # d/dt (i_b2_q) = 0
    # Bus 3
    derivatives[I_B3_D_IDX] = i_3_d + i_23_d + i_13_d - i_b3_d                                                # d/dt (i_b3_d) = 0
    derivatives[I_B3_Q_IDX] = i_3_q + i_23_q + i_13_q - i_b3_q                                                # d/dt (i_b3_q) = 0

    # Compute apparent power to return (Milano Eigenvalue Problems Eq. 1.42)
    P_terminal = v_1_d * i_1_d + v_1_q * i_1_q
    Q_terminal = v_1_q * i_1_d - v_1_d * i_1_q
    S_terminal = complex(P_terminal, Q_terminal)

    # Compute positive sequence terminal voltage to return (Milano Eigenvalue Problems Eq 1.43)
    V_terminal = (v_1_d + im * v_1_q) * exp(-im * π / 2)

    # Compute positive sequence terminal injected current to return (can check this against Milano Eigenvalue Problems Eq 1.41)
    I_terminal = conj(S_terminal / V_terminal)

    # Return all values (used for debugging – we can return fewer values after we figure out what's wrong)
    return V_terminal, S_terminal, I_terminal, i_1_d, i_1_q
end
end # module
