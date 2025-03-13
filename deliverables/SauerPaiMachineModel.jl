# Defining the module
module SauerPaiMachineModel
# Exporting the required functions and objects
export SauerPaiMachine, initialize_machine, update_machine_states!

using LinearAlgebra
using DifferentialEquations

# Defining the convention for each states
const PSI_Q = 1
const PSI_D = 2
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
    Xd::Float64
    Xq::Float64
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
    
    # Constructor with default values
    function SauerPaiMachine(;
        R = 0.002, 
        Xd = 1.79, 
        Xq = 1.71, 
        Xd_p = 0.169, 
        Xq_p = 0.228, 
        Xd_pp = 0.135,
        Xq_pp = 0.2, 
        Xl = 0.13,
        Td0_p = 4.3, 
        Tq0_p = 0.85, 
        Td0_pp = 0.032, 
        Tq0_pp = 0.05,
        γ_d1 = 0.0, 
        γ_q1 = 0.0, 
        γ_d2 = 0.0, 
        γ_q2 = 0.0,
        base_power = 100.0,
        system_base_power = 100.0,
        system_base_frequency = 60.0
    )
        return new(R, Xd, Xq, Xd_p, Xq_p, Xd_pp, Xq_pp, Xl, 
                  Td0_p, Tq0_p, Td0_pp, Tq0_pp, 
                  γ_d1, γ_q1, γ_d2, γ_q2, 
                  base_power, system_base_power, system_base_frequency)
    end
end

```
The module has two functions:
Initializing states, and updating states
This will help us in scalability
```

# Initialize machine states based on power flow solution
function initialize_machine(machine::SauerPaiMachine, V_terminal, delta, P, Q)
    # Calculate initial values for the six states based on steady-state conditions
    # This is a simplified initialization

    V_mag = abs(V_terminal)
    I_complex = conj((P + im*Q) / V_terminal)
    I_mag = abs(I_complex)
    power_factor_angle = angle(V_terminal) - angle(I_complex)
    
    # Convert terminal voltage to dq frame
    V_dq = ri_dq(delta) * [real(V_terminal); imag(V_terminal)]
    
    # Initial currents in dq frame
    I_dq = ri_dq(delta) * [real(I_complex); imag(I_complex)]
    i_d = I_dq[1]
    i_q = I_dq[2]
    
    # Initialize states - simplified for steady state
    ψd = V_dq[1] + machine.R * i_d
    ψq = V_dq[2] + machine.R * i_q
    
    eq_p = ψd + (machine.Xd - machine.Xd_p) * i_d
    ed_p = ψq - (machine.Xq - machine.Xq_p) * i_q
    
    ψd_pp = eq_p - (machine.Xd_p - machine.Xl) * i_d
    ψq_pp = -ed_p - (machine.Xq_p - machine.Xl) * i_q
    
    # Field voltage calculation
    Vf = eq_p + (machine.Xd - machine.Xd_p) * i_d
    
    # Return initial states and field voltage
    return [ψq, ψd, eq_p, ed_p, ψd_pp, ψq_pp], Vf
end

# Update machine states 
function update_machine_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    delta::Float64,
    omega::Float64,
    Vf::Float64,
    V_terminal::Complex{Float64},
    machine::SauerPaiMachine
)
    # Extract machine states
    ψq = states[PSI_Q]
    ψd = states[PSI_D]
    eq_p = states[EQ_P]
    ed_p = states[ED_P]
    ψd_pp = states[PSI_D_PP]
    ψq_pp = states[PSI_Q_PP]
    
    # System base values
    f0 = machine.system_base_frequency
    
    # Terminal voltage in dq reference frame
    V_dq = ri_dq(delta) * [real(V_terminal); imag(V_terminal)]
    

    # Calculate currents
    # Equations 15.13 in Milano's book
    i_d = (1.0 / machine.Xd_pp) * (machine.γ_d1 * eq_p - ψd + (1 - machine.γ_d1) * ψd_pp)
    i_q = (1.0 / machine.Xq_pp) * (-machine.γ_q1 * ed_p - ψq + (1 - machine.γ_q1) * ψq_pp)
    
    # Electromagnetic torque calculation (equation 15.6 in Milano's book)
    τ_e = ψd * i_q - ψq * i_d
    
    # State derivatives
    # equation 15.6 of Milano's book
    derivatives[PSI_Q] = 2.0 * π * f0 * (machine.R * i_q - omega * ψd + V_dq[2])
    derivatives[PSI_D] = 2.0 * π * f0 * (machine.R * i_d + omega * ψq + V_dq[1])
    
    # flux equations (15.13 in Milano's book)
    derivatives[EQ_P] = (1.0 / machine.Td0_p) * (
        -eq_p + Vf - (machine.Xd - machine.Xd_p) * (
            i_d - machine.γ_d2 * ψd_pp - (1 - machine.γ_d1) * i_d + machine.γ_d2 * eq_p
        )
    )
    
    # Also 15.13 in Milano's book
    derivatives[ED_P] = (1.0 / machine.Tq0_p) * (
        -ed_p + (machine.Xq - machine.Xq_p) * (
            i_q - machine.γ_q2 * ψq_pp - (1 - machine.γ_q1) * i_q - machine.γ_q2 * ed_p
        )
    )
    
    # Subtransient flux derivatives
    derivatives[PSI_D_PP] = (1.0 / machine.Td0_pp) * (eq_p - ψd_pp - (machine.Xd_p - machine.Xl) * i_d)
    derivatives[PSI_Q_PP] = (1.0 / machine.Tq0_pp) * (-ed_p - ψq_pp - (machine.Xq_p - machine.Xl) * i_q)
    
    # Calculate grid current
    I_dq = [i_d; i_q]
    I_RI = (machine.base_power / machine.system_base_power) * dq_ri(delta) * I_dq
    I_grid = Complex(I_RI[1], I_RI[2])
    
    return τ_e, I_grid
end
end # module 