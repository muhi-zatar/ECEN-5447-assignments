# Defining the module
module SauerPaiMachineModel
# Exporting the required functions and objects
export SauerPaiMachine, initialize_machine, update_machine_states!

using LinearAlgebra
using DifferentialEquations

# Defining the convention for each states
const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6

# Helper functions
# not all of them are used.

# RI to DQ transformation matrix
function ri_dq(delta)
    return [cos(delta) sin(delta); -sin(delta) cos(delta)]
end

# DQ to RI transformation matrix
function dq_ri(delta)
    return [cos(delta) -sin(delta); sin(delta) cos(delta)]
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
        γ_d1=(Xd_pp - Xl) / (Xd_p - Xl), # Removed "machine" as the object is not yet created.
        γ_q1=(Xq_pp - Xl) / (Xq_p - Xl), # Removed "machine" as the object is not yet created.
        γ_d2=(1 - γ_d1) / (Xd_p - Xl), # Removed "machine" as the object is not yet created.
        γ_q2=(1 - γ_q1) / (Xq_p - Xl), # Removed "machine" as the object is not yet created.
        H=3.148,
        D=2.0,
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

        return new(R, X_d, X_q, Xd_p, Xq_p, Xd_pp, Xq_pp, Xl,
            Td0_p, Tq0_p, Td0_pp, Tq0_pp,
            γ_d1, γ_q1, γ_d2, γ_q2, base_power,
            system_base_power, system_base_frequency, H, D)
    end
end

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
    eq_p = solution[1]
    ed_p = solution[2]
    ψd_pp = solution[3]
    ψq_pp = solution[4]
    Vf = solution[5]
    # (eq_p, ed_p, ψd_pp, ψq_pp, Vf) = A \ b

    # Use Milano Eqn. 15.11 to reconstruct the synchronous fluxes
    ψ_d = machine.R * i_q + v_q
    ψ_q = -1 * (machine.R * i_d + v_d)

    # Use Milano Eqn. 15.6 to reconstruct the electrical torque
    τ_e = ψ_d * i_q - ψ_q * i_d

    # In steady-state, τ_e = τ_m
    τ_m = τ_e

    # Prepare state vector for return
    states = zeros(Float64, 6)
    states[DELTA] = δ
    states[OMEGA] = 1.0
    states[EQ_P] = eq_p
    states[ED_P] = ed_p
    states[PSI_D_PP] = ψd_pp
    states[PSI_Q_PP] = ψq_pp

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

    # System base values
    f0 = machine.system_base_frequency

    # Extract bus voltage angle and magnitude from complex terminal voltage
    v_mag = abs(V_terminal)
    θ = angle(V_terminal)

    # Move from network to machine DQ reference frame
    # Terminal voltage in dq reference frame
    v_d = v_mag * sin(δ - θ)       # Eqn 15.4
    v_q = v_mag * cos(δ - θ)       # Eqn 15.4

    # Calculate currents (Use 15.11 and 15.15 to eliminate ψq and ψd, solve for Id and Iq)
    A = [machine.Xd_pp machine.R; machine.R -machine.Xq_pp]
    b = [-v_q + machine.γ_d1 * eq_p + (1 - machine.γ_d1) * ψd_pp; -v_d + machine.γ_q1 * ed_p - (1 - machine.γ_q1) * ψq_pp]
    currents = A \ b
    i_d = currents[1]
    i_q = currents[2]

    # We need synchronous fluxes to get electrical torque
    ψ_d = machine.R * i_q + v_q
    ψ_q = -1 * (machine.R * i_d + v_d)
    τ_e = ψ_d * i_q - ψ_q * i_d


    # State derivatives
    # shaft equations (15.5 in Milano's book)
    ω_sys = 1.0
    derivatives[DELTA] = 2.0 * π * f0 * (ω - ω_sys)

    # Speed derivative
    derivatives[OMEGA] = 1.0 / (2.0 * machine.H) * (
        τ_m - τ_e - (machine.D * (ω - 1.0))
    )


    # flux equations (15.13 in Milano's book)
    derivatives[EQ_P] = (1.0 / machine.Td0_p) * (
        -eq_p + Vf - (machine.X_d - machine.Xd_p) * (
            i_d - machine.γ_d2 * ψd_pp - (1 - machine.γ_d1) * i_d + machine.γ_d2 * eq_p
        )
    )

    # Also 15.13 in Milano's book
    derivatives[ED_P] = (1.0 / machine.Tq0_p) * (
        -ed_p + (machine.X_q - machine.Xq_p) * (
            i_q - machine.γ_q2 * ψq_pp - (1 - machine.γ_q1) * i_q - machine.γ_q2 * ed_p
        )
    )

    # Subtransient flux derivatives
    derivatives[PSI_D_PP] = (1.0 / machine.Td0_pp) * (eq_p - ψd_pp - (machine.Xd_p - machine.Xl) * i_d)
    derivatives[PSI_Q_PP] = (1.0 / machine.Tq0_pp) * (-ed_p - ψq_pp - (machine.Xq_p - machine.Xl) * i_q)

    # TODO
    # Calculate grid current
    I_dq = [i_d; i_q]
    I_RI = (machine.base_power / machine.system_base_power) * dq_ri(δ) * I_dq
    I_grid = Complex(I_RI[1], I_RI[2])

    # TODO: Convert machine DQ reference frame back to network to be good friends
    I_network_dq = [0.0; 0.0]
    V_network_dq = [0.0; 0.0]

    # TODO: After incorporating network model, return network V,I instead of machine V,I
    return V_dq, I_dq
end
end # module 