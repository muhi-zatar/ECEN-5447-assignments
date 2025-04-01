# Defining the module
module NetworkModel

# Exporing the functions
export ThreeBusNetwork, initialize_network, update_network_states!, BUS_MACHINE_MODEL, sanity_check

using LinearAlgebra

# Definining some constants
const BUS_INFINITE_BUS = 1
const BUS_MACHINE_MODEL = 2
const BUS_LOAD = 3
const NUM_STATES = 12
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
# const I_1_D_IDX = 13
# const I_1_Q_IDX = 14
# const I_3_D_IDX = 15
# const I_3_Q_IDX = 16

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
    # return Float64[
    #     sin(δ) -cos(δ)
    #     cos(δ) sin(δ)
    # ]
    return [
        0.0 -1.0
        1.0 0.0
    ]
end

function dq_ri(δ::T) where {T<:Number}
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    # return Float64[
    #     sin(δ) cos(δ)
    #     -cos(δ) sin(δ)
    # ]
    return [
        0.0 1.0
        -1.0 0.0
    ]
end

function dq_transformer(voltage, current, voltage_angle, current_angle)
    # Print initial values
    println("Voltage = $voltage")
    println("Current = $current")
    #println("Angle = $angle")

    # First, calculate with ri_dq function
    v_dq = ri_dq(voltage_angle) * voltage
    i_dq = ri_dq(current_angle) * current
    println("DQ Voltage = $v_dq")
    println("DQ Current = $i_dq")

    #####
    # Try the reverse
    v_ri = dq_ri(voltage_angle) * v_dq
    i_ri = dq_ri(current_angle) * i_dq
    println("RI Voltage = $v_ri")
    println("RI Current = $i_ri")

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

    # Reconstruct complex voltage and current at the terminal to calculate DQ
    S = complex.(P, Q)
    V_terminal = V_m .* ℯ .^ (im .* θ)
    print("\n\nInitial positive sequence voltage at the bus: $(V_terminal[2])\n\n")
    I_terminal = conj.(S ./ V_terminal)
    β = angle.(I_terminal)
    print("\n\nInitial positive sequence current at the bus: $(I_terminal[2])\n\n")

    ##### MOVE TO NETWORK DQ REFERENCE FRAME #####
    # Find DQ voltage (Milano Eigenvalue Problems Eqn 1.44)
    I_dq = zeros(Complex{Float64}, 3)
    V_dq = zeros(Complex{Float64}, 3)
    for i in 1:3
        V_RI = [real(V_terminal[i]); imag(V_terminal[i])]
        I_RI = [real(I_terminal[i]); imag(I_terminal[i])]
        V_DQ, I_DQ = dq_transformer(V_RI, I_RI, θ[i], β[i])
        V_dq[i] = complex(V_DQ[1], V_DQ[2])
        I_dq[i] = complex(I_DQ[1], I_DQ[2])
    end
    i_d = real.(I_dq)
    i_q = imag.(I_dq)
    v_d = real.(V_dq)
    v_q = imag.(V_dq)

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
    i_1_d = i_d[1]                                      # Bus 1 d-axis injected current
    i_1_q = i_q[1]                                      # Bus 1 q-axis injected current
    i_3_d = i_d[3]                                      # Bus 3 d-axis injected current
    i_3_q = i_q[3]                                      # Bus 3 q-axis injected current
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
    network.E_IB_D = v_d[BUS_INFINITE_BUS] - i_1_q * network.X_IB
    network.E_IB_Q = v_q[BUS_INFINITE_BUS] + i_1_d * network.X_IB

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

    # DEBUGGING – PRINT LINE CURRENTS
    I_12 = dq_ri(0) * [line_currents[1]; line_currents[2]]
    I_13 = dq_ri(0) * [line_currents[3]; line_currents[4]]
    I_23 = dq_ri(0) * [line_currents[5]; line_currents[6]]
    println("I_12 = $I_12")
    println("I_13 = $I_13")
    println("I_23 = $I_23")

    # Sanity check on line currents (this is just a KVL on the calculated initial values)
    res_i = [
        i_1_d - states[I_12_D_IDX] - states[I_13_D_IDX] - i_b1_d;
        i_1_q - states[I_12_Q_IDX] - states[I_13_Q_IDX] - i_b1_q;
        i_2_d + states[I_12_D_IDX] - states[I_23_D_IDX] - i_b2_d;
        i_2_q + states[I_12_Q_IDX] - states[I_23_Q_IDX] - i_b2_q;
        i_3_d + states[I_23_D_IDX] + states[I_13_D_IDX] - i_b3_d;
        i_3_q + states[I_23_Q_IDX] + states[I_13_Q_IDX] - i_b3_q
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
    # i_1_d = states[I_1_D_IDX]
    # i_1_q = states[I_1_Q_IDX]
    # i_3_d = states[I_3_D_IDX]
    # i_3_q = states[I_3_Q_IDX]
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

    # Unpack real and imaginary components of the load for convenience
    R_load = real(network.Z_L)
    X_load = imag(network.Z_L)

    # Injected currents (inner variable)
    # Bus 1 (Infinite Bus) (algebraic) –– From Milano's Eigenvalue Problems book, Eq. 1.50
    i_1_q = -1 * (network.E_IB_D - v_1_d) / network.X_IB                                  # d/dt (i_1_d) = 0
    i_1_d = (network.E_IB_Q - v_1_q) / network.X_IB                                       # d/dt (i_1_q) = 0

    # Bus 3 (Load) (inner variable) –– From Milano's Eigenvalue Problems book, Eq. 1.50
    A = [R_load -1*X_load;
        X_load R_load]
    b = [v_3_d; v_3_q]
    injected_currents = A \ b
    i_3_d = injected_currents[1]
    i_3_q = injected_currents[2]

    # Shunt currents (inner variable)
    # Bus 1
    i_b1_d = i_1_d - i_12_d - i_13_d                                                # d/dt (i_b1_d) = 0
    i_b1_q = i_1_q - i_12_q - i_13_q                                                # d/dt (i_b1_q) = 0
    # Bus 2
    i_b2_d = i_2_d + i_12_d - i_23_d                                                # d/dt (i_b2_d) = 0
    i_b2_q = i_2_q + i_12_q - i_23_q                                                # d/dt (i_b2_q) = 0
    # Bus 3
    i_b3_d = i_3_d + i_23_d + i_13_d                                                 # d/dt (i_b3_d) = 0
    i_b3_q = i_3_q + i_23_q + i_13_q                                                 # d/dt (i_b3_q) = 0

    # print the above values
    # println("i_b1_q: $i_b1_q")
    # println("i_b1_d: $i_b1_d")
    # println("i_b2_q: $i_b2_q")
    # println("i_b2_d: $i_b2_d")
    # println("i_b3_q: $i_b3_q")
    # println("i_b3_d: $i_b3_d")
    # exit()
    # Define derivatives of all states (differential and algebraic)
    # Line currents (differential)
    derivatives[I_12_D_IDX] = (v_1_d - v_2_d - network.R_12 * i_12_d + network.X_12 * i_12_q) / network.X_12                 # d/dt (i_12_d) != 0
    derivatives[I_12_Q_IDX] = (v_1_q - v_2_q - network.R_12 * i_12_q - network.X_12 * i_12_d) / network.X_12                 # d/dt (i_12_q) != 0
    derivatives[I_13_D_IDX] = (v_1_d - v_3_d - network.R_13 * i_13_d + network.X_13 * i_13_q) / network.X_13                 # d/dt (i_13_d) != 0
    derivatives[I_13_Q_IDX] = (v_1_q - v_3_q - network.R_13 * i_13_q - network.X_13 * i_13_d) / network.X_13                 # d/dt (i_13_q) != 0
    derivatives[I_23_D_IDX] = (v_2_d - v_3_d - network.R_23 * i_23_d + network.X_23 * i_23_q) / network.X_23                 # d/dt (i_23_d) != 0
    derivatives[I_23_Q_IDX] = (v_2_q - v_3_q - network.R_23 * i_23_q - network.X_23 * i_23_d) / network.X_23                 # d/dt (i_23_q) != 0

    # derivatives[I_12_D_IDX] = 0
    # derivatives[I_12_Q_IDX] = 0
    # derivatives[I_13_D_IDX] = 0
    # derivatives[I_13_Q_IDX] = 0
    # derivatives[I_23_D_IDX] = 0
    # derivatives[I_23_Q_IDX] = 0
    # Bus voltages (differential)
    derivatives[V_1_D_IDX] = (i_b1_d + network.B_1 * v_1_q) / network.B_1                                                     # d/dt (v_1_d) != 0
    derivatives[V_1_Q_IDX] = (i_b1_q - network.B_1 * v_1_d) / network.B_1                                                     # d/dt (v_1_q) != 0
    derivatives[V_2_D_IDX] = (i_b2_d + network.B_2 * v_2_q) / network.B_2                                                     # d/dt (v_2_d) != 0
    derivatives[V_2_Q_IDX] = (i_b2_q - network.B_2 * v_2_d) / network.B_2                                                     # d/dt (v_2_q) != 0
    derivatives[V_3_D_IDX] = (i_b3_d + network.B_3 * v_3_q) / network.B_3                                                     # d/dt (v_3_d) != 0
    derivatives[V_3_Q_IDX] = (i_b3_q - network.B_3 * v_3_d) / network.B_3                                                     # d/dt (v_3_q) != 0

    # put the above derivates to zero
    # derivatives[V_1_D_IDX] = 0
    # derivatives[V_1_Q_IDX] = 0
    # derivatives[V_2_D_IDX] = 0
    # derivatives[V_2_Q_IDX] = 0
    # derivatives[V_3_D_IDX] = 0
    # derivatives[V_3_Q_IDX] = 0
    # Compute apparent power to return (Milano Eigenvalue Problems Eq. 1.42)
    P_terminal = v_2_d * i_2_d + v_2_q * i_2_q
    Q_terminal = v_2_q * i_2_d - v_2_d * i_2_q
    S_terminal = complex(P_terminal, Q_terminal)
    #sanity_check(S_terminal, S_bus, "Machine vs. Network Power", false)

    # Compute positive sequence terminal voltage to return (Milano Eigenvalue Problems Eq 1.43)
    V_terminal = (v_2_d + im * v_2_q) * exp(-im * π / 2)

    # Compute positive sequence terminal injected current to return (can check this against Milano Eigenvalue Problems Eq 1.41)
    I_terminal = conj(S_terminal / V_terminal)
    I_terminal_test = (i_2_d + im * i_2_q) * exp(-im * π / 2)
    #sanity_check(I_terminal_test, I_terminal, "Bus 2 Current Injection", false)

    # Return all values (used for debugging – we can return fewer values after we figure out what's wrong)
    return V_terminal, S_terminal, I_terminal, i_2_d, i_2_q
end
end # module