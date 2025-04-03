module GovernorComponents

export GasTG, initialize_gov_states, update_gov_states!
export FV_IDX, FF_IDX, ET_IDX

# Governor state indices
const FV_IDX = 1    # Fuel valve
const FF_IDX = 2    # Fuel flow
const ET_IDX = 3    # Exhaust temp

"""
    GasTG

Gas turbine governor model with parameters.
"""
mutable struct GasTG
    # Governor parameters
    R::Float64       # Droop (pu)
    T1::Float64      # Governor time constant (s)
    T2::Float64      # Combustor time constant (s)
    T3::Float64      # Load limiter time constant (s)
    D_turb::Float64  # Turbine damping factor (pu)
    AT::Float64      # Ambient temperature load limit (pu)
    KT::Float64      # Temperature control loop gain
    V_min::Float64   # Minimum valve position (pu)
    V_max::Float64   # Maximum valve position (pu)
    P_ref::Float64   # Reference power set point (pu)

    # Constructor with default values
    function GasTG(;
        R=0.05,
        T1=0.2,
        T2=0.2,
        T3=2.0,
        D_turb=0.0,
        AT=1.0,
        KT=2.5,
        V_min=0.01,
        V_max=1.1,
        P_ref=1.0
    )
        return new(R, T1, T2, T3, D_turb, AT, KT, V_min, V_max, P_ref)
    end
end

"""
    initialize_gov_states(gov::GasTG, τ_m_init::Float64, ω_init::Float64)

Initialize governor states based on the steady-state rotor angular velocity 
and the steady-state mechanical torque.
"""
function initialize_gov_states(gov::GasTG, τ_m_init::Float64, ω_init::Float64)
    # This function initializes the governor states based on steady-state conditions
    states = zeros(Float64, 3)

    # Initialize all states to the initial torque value (steady-state)
    states[FV_IDX] = τ_m_init  # Fuel valve position
    states[FF_IDX] = τ_m_init  # Fuel flow
    states[ET_IDX] = τ_m_init  # Exhaust temperature

    return states
end

"""
    update_gov_states!(states, derivatives, omega, gov)

Update governor states using state space formulation.
Returns the mechanical torque.
"""
function update_gov_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    omega::Float64,
    gov::GasTG
)
    # Extract state variables
    fv = states[FV_IDX]  # Fuel valve position
    ff = states[FF_IDX]  # Fuel flow
    et = states[ET_IDX]  # Exhaust temperature

    # Calculate speed error (delta omega)
    # Δω = ω - ωref where ωref = 1.0 (per unit)
    delta_omega = omega - 1.0

    # Calculate inverse droop (1/R) with protection against division by zero
    inv_R = gov.R < eps() ? 0.0 : 1.0 / gov.R

    # Calculate inputs for the two control paths
    # Power reference path: Pref - (1/R)Δω
    power_reference_input = gov.P_ref - inv_R * delta_omega

    # Temperature limit path: AT + KT(AT - ET)
    temperature_limit_input = gov.AT + gov.KT * (gov.AT - et)

    # Determine active path using MIN function from the diagram
    governor_input = min(power_reference_input, temperature_limit_input)

    # Apply valve position limits if needed
    governor_input = max(gov.V_min, min(gov.V_max, governor_input))

    # State equations for each component
    
    # Fuel flow derivative (first integrator with time constant T₁)
    dff_dt = (governor_input - ff) / gov.T1

    # Fuel valve derivative (second integrator with time constant T₂)
    dfv_dt = (ff - fv) / gov.T2

    # Exhaust temperature derivative (third integrator with time constant T₃)
    det_dt = (fv - et) / gov.T3

    # Update derivatives vector
    derivatives[FV_IDX] = dfv_dt
    derivatives[FF_IDX] = dff_dt
    derivatives[ET_IDX] = det_dt

    # Calculate mechanical power output
    # Pm = x₁ - D_turb·Δω (fuel valve position minus damping term)
    Pm = fv - gov.D_turb * delta_omega

    # Convert from power to torque for the mechanical system (P = ω * τ, so τ = P/ω)
    τm = Pm / omega

    return τm
end

end # module