# Defining the module
module SauerPaiMachineModel
# Exporting the required functions and objects
export SauerPaiMachine, initialize_machine, update_machine_states!

using LinearAlgebra
using DifferentialEquations
using NLsolve
# Defining the convention for each states
const NUM_STATES = 8
const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6
const PSI_D = 7
const PSI_Q = 8

const STRICT_NLSOLVE_F_TOLERANCE = 1e-8

# Helper functions
# not all of them are used.

# RI to DQ transformation matrix
function ri_dq_machine(δ::T) where {T<:Number}
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return Float64[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

function dq_ri_machine(δ::T) where {T<:Number}
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return Float64[
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end


# Machine structure for parameters
mutable struct SauerPaiMachine
    # Machine parameters
    R::Float64
    X_d::Float64
    X_q::Float64
    Xd_p::Float64
    Xq_p::Float64
    Xd_pp::Float64
    Xq_pp::Float64
    Xl::Float64
    Td0_p::Float64
    Tq0_p::Float64
    Td0_pp::Float64
    Tq0_pp::Float64
    γ_d1::Float64
    γ_q1::Float64
    γ_d2::Float64
    γ_q2::Float64
    base_power::Float64
    system_base_power::Float64
    system_base_frequency::Float64
    H::Float64 # Rearragned for tideness
    D::Float64 # REarranged for tindeness
    M::AbstractArray{Float64}           # To be compatible with mass-matrix formulation (use 1s and 0s)

    # Constructor with default values
    function SauerPaiMachine(;
        R=0.002,
        X_d=1.79,
        X_q=1.71,
        Xd_p=0.169,
        Xq_p=0.228,
        Xd_pp=0.135,
        Xq_pp=0.2,
        Xl=0.13,
        Td0_p=4.3,
        Tq0_p=0.85,
        Td0_pp=0.032,
        Tq0_pp=0.05,
        γ_d1=0.1282051282051283,
        γ_q1=0.7142857142857143,
        γ_d2=22.35371466140696,
        γ_q2=2.9154518950437316,
        H=3.148,
        # H = 10000.0,
        D=2.0,
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0,
        M=ones(Float64, NUM_STATES)
    )

        return new(R, X_d, X_q, Xd_p, Xq_p, Xd_pp, Xq_pp, Xl,
            Td0_p, Tq0_p, Td0_pp, Tq0_pp,
            γ_d1, γ_q1, γ_d2, γ_q2, base_power,
            system_base_power, system_base_frequency, H, D, M)
    end
end

get_R(machine::SauerPaiMachine) = machine.R
get_Xd(machine::SauerPaiMachine) = machine.X_d
get_Xq(machine::SauerPaiMachine) = machine.X_q
get_Xd_p(machine::SauerPaiMachine) = machine.Xd_p
get_Xq_p(machine::SauerPaiMachine) = machine.Xq_p
get_Xd_pp(machine::SauerPaiMachine) = machine.Xd_pp
get_Xq_pp(machine::SauerPaiMachine) = machine.Xq_pp
get_Xl(machine::SauerPaiMachine) = machine.Xl
get_Td0_p(machine::SauerPaiMachine) = machine.Td0_p
get_Tq0_p(machine::SauerPaiMachine) = machine.Tq0_p
get_Td0_pp(machine::SauerPaiMachine) = machine.Td0_pp
get_Tq0_pp(machine::SauerPaiMachine) = machine.Tq0_pp
get_γ_d1(machine::SauerPaiMachine) = machine.γ_d1
get_γ_q1(machine::SauerPaiMachine) = machine.γ_q1
get_γ_d2(machine::SauerPaiMachine) = machine.γ_d2
get_γ_q2(machine::SauerPaiMachine) = machine.γ_q2

```
The module has two functions:
Initializing states, and updating states
This will help us in scalability
```

# Initialize machine states based on power flow solution
function initialize_machine(machine::SauerPaiMachine, V_terminal, delta, P, Q)
    ```
    This function will initialize the states of the machine at steady-state.
    It will also initialize the field voltage, to be used by initialize_avr, and the mechanical
    torque, to be used by initialize_gov.
    ```
    # Calculate initial values for the six states based on steady-state conditions
    # This is a simplified initialization
    R = get_R(machine)
    Xd = get_Xd(machine)
    Xq = get_Xq(machine)
    Xd_p = get_Xd_p(machine)
    Xq_p = get_Xq_p(machine)
    Xd_pp = get_Xd_pp(machine)
    Xq_pp = get_Xq_pp(machine)
    Xl = get_Xl(machine)
    Td0_p = get_Td0_p(machine)
    γ_d1 = get_γ_d1(machine)
    γ_q1 = get_γ_q1(machine)
    γ_d2 = get_γ_d2(machine)
    γ_q2 = get_γ_q2(machine)

    V_mag = abs(V_terminal)
    I_complex = conj((P + im * Q) / V_terminal)
    I_mag = abs(I_complex)
    power_factor_angle = angle(V_terminal) - angle(I_complex)

    # Calculate initial rotor angle
    δ = angle(V_terminal + complex(machine.R, machine.X_q) * I_complex)             # Eqn. 9.11

    # Convert terminal voltage to dq frame
    vdq = V_terminal * ℯ^(-1 * im * (δ - π / 2))
    v_d = real(vdq)
    v_q = imag(vdq)

    # Current
    idq = I_complex * ℯ^(-1 * im * (δ - π / 2)) # Fixed from δ to δ0
    i_d = real(idq)
    i_q = imag(idq)

    # Reconstruct power and print for sanity check:
    p_bus = v_d * i_d + v_q * i_q
    q_bus = v_q * i_d - v_d * i_q
    println("Initial S_bus=$(complex(p_bus, q_bus))")       # This is correct

    # Initialize states - simplified for steady state
    A = [-machine.γ_q1 0 0 (1-machine.γ_q1) 0;
        0 machine.γ_d1 (1-machine.γ_d1) 0 0;
        0 1 -1 0 0;
        -1 0 0 -1 0;
        0 -1 0 0 1]
    b = [machine.Xq_pp * i_q - machine.R * i_d - v_d;
        machine.R * i_q + machine.Xd_pp * i_d + v_q;
        (machine.Xd_p - machine.Xl) * i_d;
        (machine.Xq_p - machine.Xl) * i_q;
        (machine.X_d - machine.Xd_p) * i_d]

    # Tideness and easier tracability
    solution = A \ b
    ed_p = solution[1]
    eq_p = solution[2]
    ψd_pp = solution[3]
    ψq_pp = solution[4]
    Vf = solution[5]
    # (eq_p, ed_p, ψd_pp, ψq_pp, Vf) = A \ b

    # Use Milano Eqn. 15.11 to reconstruct the synchronous fluxes
    ψ_d = machine.R * i_q + v_q
    ψ_q = -1 * (machine.R * i_d + v_d)

    # Use Milano Eqn. 15.6 to reconstruct the electrical torque
    τ_e = ψ_d * i_q - ψ_q * i_d
    ω0 = 1.0
    # In steady-state, τ_e = τ_m
    τ_m = τ_e
    V_R = real(V_terminal)
    V_I = imag(V_terminal)
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        ψq = x[4]
        ψd = x[5]
        eq_p = x[6]
        ed_p = x[7]
        ψd_pp = x[8]
        ψq_pp = x[9]

        # Milano
        V_dq = ri_dq_machine(δ) * [V_R; V_I]
        i_d = (1.0 / Xd_pp) * (γ_d1 * eq_p - ψd + (1 - γ_d1) * ψd_pp)      #15.15
        i_q = (1.0 / Xq_pp) * (-γ_q1 * ed_p - ψq + (1 - γ_q1) * ψq_pp)     #15.15
        τ_e = ψd * i_q - ψq * i_d               #15.6
        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P - (V_dq[1] * i_d + V_dq[2] * i_q) #Output Power (15.2)
        out[3] = Q - (V_dq[2] * i_d - V_dq[1] * i_q) #Output Reactive Power (15.3)
        out[4] = R * i_q - ω0 * ψd + V_dq[2]                                    #15.9 ψq
        out[5] = R * i_d + ω0 * ψq + V_dq[1]                                    #15.9 ψd
        out[6] =
            -eq_p - (Xd - Xd_p) * (i_d - γ_d2 * ψd_pp - (1 - γ_d1) * i_d + γ_d2 * eq_p) +
            Vf0    #15.13
        out[7] = -ed_p + (Xq - Xq_p) * (i_q - γ_q2 * ψq_pp - (1 - γ_q1) * i_q - γ_d2 * ed_p)          #15.13
        out[8] = -ψd_pp + eq_p - (Xd_p - Xl) * i_d      #15.13
        out[9] = -ψq_pp - ed_p - (Xq_p - Xl) * i_q #15.13
    end

    x0 = [δ, τ_m, Vf, ψ_q, ψ_d, eq_p, ed_p, ψd_pp, ψq_pp]
    sol = NLsolve.nlsolve(f!, x0; ftol=STRICT_NLSOLVE_F_TOLERANCE)

    if !NLsolve.converged(sol)
        @warn("Initialization in machine failed!")
    else
        # Get the solution
        sol_x0 = sol.zero

        states = zeros(Float64, 8)
        τ_m = sol_x0[2]
        Vf = sol_x0[3]
        states[DELTA] = δ
        states[OMEGA] = 1.0
        states[EQ_P] = eq_p
        states[ED_P] = ed_p
        states[PSI_D_PP] = ψd_pp
        states[PSI_Q_PP] = ψq_pp
        states[PSI_D] = ψ_d
        states[PSI_Q] = ψ_q
    end

    # Return initial states and field voltage
    return states, Vf, τ_m
end

# Update machine states 
function update_machine_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    V_terminal::Complex{Float64},
    Vf::Float64,
    τ_m::Float64,
    machine::SauerPaiMachine
)
    # Extract machine states
    δ = states[DELTA]
    ω = states[OMEGA]
    eq_p = states[EQ_P]
    ed_p = states[ED_P]
    ψd_pp = states[PSI_D_PP]
    ψq_pp = states[PSI_Q_PP]
    ψ_d = states[PSI_D]
    ψ_q = states[PSI_Q]
    # System base values
    f0 = machine.system_base_frequency


    V_dq = ri_dq_machine(δ) * [real(V_terminal); imag(V_terminal)]
    v_d = V_dq[1]       # Eqn 15.4
    v_q = V_dq[2]       # Eqn 15.4
    V_mag = abs(complex(v_d, v_q)) # Debugging

    i_d = (-ψ_d + machine.γ_d1 * eq_p + (1 - machine.γ_d1) * ψd_pp) / machine.Xd_pp
    i_q = (-ψ_q - machine.γ_q1 * ed_p + (1 - machine.γ_q1) * ψq_pp) / machine.Xq_pp
    # Calculate electrical torque
    τ_e = ψ_d * i_q - ψ_q * i_d
    I_mag = abs(complex(i_d, i_q)) # Debugging

    # State derivatives
    # shaft equations (15.5 in Milano's book)
    # Fixing the system frequency to machine frequency
    ω_sys = ω
    derivatives[DELTA] = 2.0 * π * f0 * (ω - ω_sys)

    # Speed derivative
    derivatives[OMEGA] = 1.0 / (2.0 * machine.H) * (
        τ_m - τ_e - (machine.D * (ω - 1.0))
    )

    # flux equations (15.9 in Milano's book)
    # These equations here are differential states
    derivatives[PSI_D] = (machine.R * i_d + ω * ψ_q + v_d) * (2.0 * π * f0)
    derivatives[PSI_Q] = (machine.R * i_q - ω * ψ_d + v_q) * (2.0 * π * f0)
    # derivatives[PSI_D] = 0.0
    # derivatives[PSI_Q] = 0.0
    # flux equations (15.13 in Milano's book)
    derivatives[EQ_P] = (1.0 / machine.Td0_p) * (
        -eq_p + Vf - (machine.X_d - machine.Xd_p) * (
            i_d - machine.γ_d2 * ψd_pp - (1 - machine.γ_d1) * i_d + machine.γ_d2 * eq_p
        )
    )

    # Also 15.13 in Milano's book
    derivatives[ED_P] = (1.0 / machine.Tq0_p) * (
        -ed_p + (machine.X_q - machine.Xq_p) * (
            i_q - machine.γ_q2 * ψq_pp - (1 - machine.γ_q1) * i_q - machine.γ_d2 * ed_p
        )
    )

    # Subtransient flux derivatives
    derivatives[PSI_D_PP] = (1.0 / machine.Td0_pp) * (eq_p - ψd_pp - (machine.Xd_p - machine.Xl) * i_d)
    derivatives[PSI_Q_PP] = (1.0 / machine.Tq0_pp) * (-ed_p - ψq_pp - (machine.Xq_p - machine.Xl) * i_q)

    # Calculate power at the bus to return
    p_bus = v_d * i_d + v_q * i_q       # Milano Eq. 15.2
    q_bus = v_q * i_d - v_d * i_q       # Milano Eq. 15.3
    S_bus = complex(p_bus, q_bus)

    # Calculate current in positive sequence
    I_RI = conj(S_bus / V_terminal)

    return I_RI, S_bus, ω, V_mag, I_mag
end
end # module 
