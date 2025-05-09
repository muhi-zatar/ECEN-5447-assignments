# Defining the module
module NetworkModel_Multimachine

# Exporing the functions
export ThreeBusNetwork, initialize_network, update_network_states!, BUS_MACHINE_MODEL, sanity_check

using LinearAlgebra

# Definining some constants
const BUS_MACHINE_MODEL = 1
const BUS_CONVERTER_MODEL = 2
const BUS_LOAD = 3
const NUM_STATES = 20
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
const I_3_D_IDX = 13
const I_3_Q_IDX = 14
const I_B1_D_IDX = 15
const I_B1_Q_IDX = 16
const I_B2_D_IDX = 17
const I_B2_Q_IDX = 18
const I_B3_D_IDX = 19
const I_B3_Q_IDX = 20

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
    machine_angle::Float64


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
        M=zeros(Float64, NUM_STATES, NUM_STATES),  # To be filled in during initialization
        machine_angle=0.0
    )
        return new(R_12, X_12, B_1, R_13, X_13, B_3, R_23, X_23, B_2, X_IB, E_IB_D, E_IB_Q, Z_L, M, machine_angle)
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
function initialize_network(network::ThreeBusNetwork, V_m::Vector{Float64},
    θ::Vector{Float64}, P::Vector{Float64},
    Q::Vector{Float64}, δ_machine::Float64=0.0)
    # Prepare state vector for return
    states = zeros(Float64, NUM_STATES)

    # Reconstruct complex voltage and current at the terminal to calculate DQ
    S = complex.(P, Q)
    V_terminal = V_m .* ℯ .^ (im .* θ)
    print("\n\nInitial positive sequence voltage at the bus: $(V_terminal[2])\n\n")
    I_terminal = conj.(S ./ V_terminal)
    print("\n\nInitial positive sequence current at the bus: $(I_terminal[2])\n\n")

    ##### MOVE TO MACHINE DQ REFERENCE FRAME INSTEAD OF NETWORK FRAME #####
    # Find DQ voltage using the machine angle as reference
    I_dq = zeros(Complex{Float64}, 3)
    V_dq = zeros(Complex{Float64}, 3)
    for i in 1:3
        V_RI = [real(V_terminal[i]); imag(V_terminal[i])]
        I_RI = [real(I_terminal[i]); imag(I_terminal[i])]
        # Use machine angle as reference frame instead of 0.0
        V_DQ, I_DQ = dq_transformer(V_RI, I_RI, δ_machine, δ_machine)
        V_dq[i] = complex(V_DQ[1], V_DQ[2])
        I_dq[i] = complex(I_DQ[1], I_DQ[2])
    end
    i_d = real.(I_dq)
    i_q = imag.(I_dq)
    v_d = real.(V_dq)
    v_q = imag.(V_dq)

    # Rest of the function remains largely unchanged...
    # Bus 1 (machine model) injection current in machine reference frame
    i_1_d = i_d[BUS_MACHINE_MODEL]
    i_1_q = i_q[BUS_MACHINE_MODEL]

    # Bus 2 (inverter model) injection current in machine reference frame
    i_2_d = i_d[BUS_CONVERTER_MODEL]
    i_2_q = i_q[BUS_CONVERTER_MODEL]

    # Sanity check
    P_test = (v_d .* i_d) .+ (v_q .* i_q)
    Q_test = (v_q .* i_d) .- (v_d .* i_q)
    S_test = complex.(P_test, Q_test)
    sanity_check(S_test, S, "Power")

    # Find load impedance - impedance is frame-invariant, so no change needed
    Z_dq = (abs.(V_dq) .^ 2) ./ conj.(S)
    network.Z_L = Z_dq[BUS_LOAD]             # The load is at the third bus


    # Define state vector elements (Assume steady-state)
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
    I_12 = dq_ri(δ_machine) * [line_currents[1]; line_currents[2]]
    I_13 = dq_ri(δ_machine) * [line_currents[3]; line_currents[4]]
    I_23 = dq_ri(δ_machine) * [line_currents[5]; line_currents[6]]
    println("I_12 = $I_12")
    println("I_13 = $I_13")
    println("I_23 = $I_23")

    # Sanity check on line currents (this is just a KVL on the calculated initial values)
    res_i = [
        i_1_d - states[I_12_D_IDX] - states[I_13_D_IDX] - states[I_B1_D_IDX];
        i_1_q - states[I_12_Q_IDX] - states[I_13_Q_IDX] - states[I_B1_Q_IDX];
        i_2_d + states[I_12_D_IDX] - states[I_23_D_IDX] - states[I_B2_D_IDX];
        i_2_q + states[I_12_Q_IDX] - states[I_23_Q_IDX] - states[I_B2_Q_IDX];
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
    M_diagonal[I_3_D_IDX] = 0.0
    M_diagonal[I_3_Q_IDX] = 0.0
    M_diagonal[I_B1_D_IDX:I_B3_Q_IDX] .= 0.0
    network.M = M_diagonal
    network.machine_angle = δ_machine

    # Return initial states
    return states, i_1_d, i_1_q, i_2_d, i_2_q
end

function update_network_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    S_machine_bus::Complex{Float64},
    S_converter_bus::Complex{Float64},
    network::ThreeBusNetwork,
    δ_machine::Float64,  # Add machine angle as parameter
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
    i_3_d = states[I_3_D_IDX]
    i_3_q = states[I_3_Q_IDX]
    i_b1_d = states[I_B1_D_IDX]
    i_b1_q = states[I_B1_Q_IDX]
    i_b2_d = states[I_B2_D_IDX]
    i_b2_q = states[I_B2_Q_IDX]
    i_b3_d = states[I_B3_D_IDX]
    i_b3_q = states[I_B3_Q_IDX]

    # Calculate machine DQ current from bus power and voltage
    V1_dq = complex(v_1_d, v_1_q)
    I1_dq = conj(S_machine_bus ./ V1_dq)
    i_1_d = real(I1_dq)
    i_1_q = imag(I1_dq)

    # Calculate converter DQ current from bus power and voltage
    V2_dq = complex(v_2_d, v_2_q)
    I2_dq = conj(S_converter_bus ./ V2_dq)
    i_2_d = real(I2_dq)
    i_2_q = imag(I2_dq)

    # Unpack real and imaginary components of the load for convenience
    R_load = real(network.Z_L)
    X_load = imag(network.Z_L)

    # Define derivatives of all states (differential and algebraic)
    # Line currents (differential)
    derivatives[I_12_D_IDX] = (v_1_d - v_2_d - network.R_12 * i_12_d + network.X_12 * i_12_q)
    derivatives[I_12_Q_IDX] = (v_1_q - v_2_q - network.R_12 * i_12_q - network.X_12 * i_12_d)
    derivatives[I_13_D_IDX] = (v_1_d - v_3_d - network.R_13 * i_13_d + network.X_13 * i_13_q)
    derivatives[I_13_Q_IDX] = (v_1_q - v_3_q - network.R_13 * i_13_q - network.X_13 * i_13_d)
    derivatives[I_23_D_IDX] = (v_2_d - v_3_d - network.R_23 * i_23_d + network.X_23 * i_23_q)
    derivatives[I_23_Q_IDX] = (v_2_q - v_3_q - network.R_23 * i_23_q - network.X_23 * i_23_d)

    # Bus voltages (differential)
    derivatives[V_1_D_IDX] = (i_b1_d + network.B_1 * v_1_q)
    derivatives[V_1_Q_IDX] = (i_b1_q - network.B_1 * v_1_d)
    derivatives[V_2_D_IDX] = (i_b2_d + network.B_2 * v_2_q)
    derivatives[V_2_Q_IDX] = (i_b2_q - network.B_2 * v_2_d)
    derivatives[V_3_D_IDX] = (i_b3_d + network.B_3 * v_3_q)
    derivatives[V_3_Q_IDX] = (i_b3_q - network.B_3 * v_3_d)

    # Current injections
    # Bus 3 (Load) (algebraic) –– From Milano's Eigenvalue Problems book, Eq. 1.48
    derivatives[I_3_D_IDX] = (v_3_d - R_load * i_3_d + X_load * i_3_q)                                          # d/dt (i_3_d) = 0
    derivatives[I_3_Q_IDX] = (v_3_q - R_load * i_3_q - X_load * i_3_d)                                          # d/dt (i_3_q) = 0

    # Shunt currents (algebraic)
    # Bus 1
    derivatives[I_B1_D_IDX] = i_1_d - i_12_d - i_13_d - i_b1_d                                                # d/dt (i_b1_d) = 0
    derivatives[I_B1_Q_IDX] = i_1_q - i_12_q - i_13_q - i_b1_q                                                # d/dt (i_b1_q) = 0
    # Debugging
    # println("I_B1_D = $i_b1_d")
    # println("I_B1_Q = $i_b1_q")
    # println("I_12_D = $i_12_d")
    # println("I_12_Q = $i_12_q")
    # println("I_13_D = $i_13_d")
    # println("I_13_Q = $i_13_q")
    # println("I_1_D = $i_1_d")
    # println("I_1_Q = $i_1_q")
    # Bus 2
    derivatives[I_B2_D_IDX] = i_2_d + i_12_d - i_23_d - i_b2_d                                                # d/dt (i_b2_d) = 0
    derivatives[I_B2_Q_IDX] = i_2_q + i_12_q - i_23_q - i_b2_q                                                # d/dt (i_b2_q) = 0
    # Bus 3
    derivatives[I_B3_D_IDX] = i_3_d + i_23_d + i_13_d - i_b3_d                                                # d/dt (i_b3_d) = 0
    derivatives[I_B3_Q_IDX] = i_3_q + i_23_q + i_13_q - i_b3_q                                                # d/dt (i_b3_q) = 0

    # Compute apparent power to return 
    P_terminal = v_1_d * i_1_d + v_1_q * i_1_q
    Q_terminal = v_1_q * i_1_d - v_1_d * i_1_q
    S_terminal = complex(P_terminal, Q_terminal)

    # For terminal voltage/current, we need to transform back to system reference frame
    # for interfacing with other components that expect system frame values
    V_terminal_machine_dq = v_1_d + im * v_1_q

    # Convert from machine DQ to system reference
    V_terminal_RI = dq_ri(δ_machine) * [v_1_d; v_1_q]
    V_terminal = V_terminal_RI[1] + im * V_terminal_RI[2]

    # Current also needs transformation
    I_terminal = conj(S_terminal / V_terminal)

    return V_terminal, S_terminal, I_terminal, i_2_d, i_2_q
end
end # module
