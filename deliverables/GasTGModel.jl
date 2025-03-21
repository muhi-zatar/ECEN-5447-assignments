# Defining the module
module GasTGModel

# Exporting the needed functions
export GasTG, initialize_gov, update_gov_states!

# Governor state indices
const FV_IDX = 1                # Fuel value
const FF_IDX = 2                # Fuel flow
const ET_IDX = 3                # Exhaust temp

# Gas turbine governor model structure
mutable struct GasTG
    # Governor parameters
    R::Float64       # in pu
    T1::Float64      # Governor time constant (s)
    T2::Float64      # Combustor time constant (s)
    T3::Float64      # Load limiter time constant (s)
    D_turb::Float64  # Turbine damping factor (pu) it is set to zero in the project description, but needed for equations
    AT::Float64      # Ambient temperature load limit (pu)
    KT::Float64      # Temperature control loop gain
    V_min::Float64   # Minimum valve position (pu)
    V_max::Float64   # Maximum valve position (pu)
    P_ref::Float64   # Reference power set point (pu)

    # Constructor with default values
    function GasTG(;
        R=0.05,
        T1=0.4,
        T2=0.1,
        T3=3.0,
        D_turb=0.0,
        AT=1.0,
        KT=2.0,
        V_min=0.0,
        V_max=1.0,
        P_ref=1.0
    )
        return new(R, T1, T2, T3, D_turb, AT, KT, V_min, V_max, P_ref)
    end
end

# Initialize governor states
function initialize_gov(gov::GasTG, τ_m_init::Float64, ω_init::Float64)
    # This function will initialize the governor states based on the steady-state rotor angular
    # velocity and the steady-state mechanical torque
    
    # Extract governor states
    states = zeros(Float64, 3)

    # Initialize all states to the initial torque value
    # This assumes steady-state initial conditions
    states[FV_IDX] = τ_m_init
    states[FF_IDX] = τ_m_init
    states[ET_IDX] = τ_m_init

    # Return initialized states
    return states
end

# Update governor states using state space formulation
function update_gov_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    omega::Float64,
    gov::GasTG
)
    # Extract state variables from the state vector:
    # x₁ = Fuel Valve (FV) - Output of the second integrator with time constant T₂
    fv = states[FV_IDX]  
    
    # x₂ = Fuel Flow (FF) - Output of the first integrator with time constant T₁
    ff = states[FF_IDX]  
    
    # x₃ = Exhaust Temperature (ET) - Output of the third integrator with time constant T₃
    et = states[ET_IDX]
    
    # Calculate speed error (delta omega)
    # Δω = ω - ωref where ωref = 1.0 (per unit)
    delta_omega = omega - 1.0
    
    # Calculate inverse droop (1/R)y
    inv_R = gov.R < eps() ? 0.0 : 1.0 / gov.R
    
    # Calculate inputs for the two paths
    # Power reference path: Pref - (1/R)Δω
    power_reference_input = gov.P_ref - inv_R * delta_omega
    
    # Temperature limit path: AT + KT(AT - ET)
    temperature_limit_input = gov.AT + gov.KT * (gov.AT - et)
    
    # Determine active path using MIN function from the diagram
    governor_input = min(power_reference_input, temperature_limit_input)
    
    # State equations expressed in standard form

    # dFF/dt equation (Fuel Flow - state x2)
    # ẋ₂ = (governor_input - x₂)/T₁
    # where governor_input = min(Pref - (1/R)Δω, AT + KT(AT - x₃))
    dff_dt = (governor_input - ff) / gov.T1
    
    # dFV/dt equation (Fuel Valve - state x1)
    # ẋ₁ = (x₂ - x₁)/T₂
    dfv_dt = (ff - fv) / gov.T2
    
    # dET/dt equation (Exhaust Temperature - state x3)
    # ẋ₃ = (x₁ - x₃)/T₃
    det_dt = (fv - et) / gov.T3
    
    # Update derivatives vector
    derivatives[FV_IDX] = dfv_dt
    derivatives[FF_IDX] = dff_dt
    derivatives[ET_IDX] = det_dt
    
    # Calculate mechanical power output
    # Pm = x₁ - D_turb·Δω
    # The mechanical power is the fuel valve position minus the damping term
    Pm = fv - gov.D_turb * delta_omega
    
    # Calculate mechanical torque (P = ω * τ, so τ = P/ω)
    # Convert from power to torque for the mechanical system
    τm = Pm / omega
    
    return τm
end

end # module