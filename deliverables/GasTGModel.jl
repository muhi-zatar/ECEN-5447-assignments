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


#The below equations are from the diagram in PowerWorld following link:
# https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20GAST_PTI%20and%20GASTD.htm



# Helper functions for the governor blocks
# When checked against Jose's implementation in the julia package, they are not the same.
# Bunch of helper functions for the AVR blocks
function low_pass(input, state, gain, time_constant)
    if time_constant > 0.0
        output = state
        derivative = (gain * input - state) / time_constant
        return derivative
    else
        # For zero time constant, output follows input
        return 0.0
    end
end

function low_pass_nonwindup(input, state, gain, time_constant, min_val, max_val)
    if time_constant > 0.0
        # Limit the state value
        limited_state = clamp(state, min_val, max_val)
        derivative = (gain * input - limited_state) / time_constant

        # Anti-windup logic
        if (limited_state >= max_val && derivative > 0.0) ||
           (limited_state <= min_val && derivative < 0.0)
            derivative = 0.0
        end

        return limited_state, derivative
    else
        # For zero time constant, output = gain * input (algebraic)
        output = clamp(gain * input, min_val, max_val)
        return output, 0.0
    end
end

function clamp(value, min_val, max_val)
    return max(min(value, max_val), min_val)
end
# END OF HELPER FUNCTIONS


# The module has two functions:
# Initializing states, and updating states
# This will help us in scalability

# Initialize governor states
function initialize_gov(gov::GasTG, τ_m_init::Float64, ω_init::Float64)
    
    # This function will initialize the governor states based on the steady-state rotor angular
    # velocity and the steady-state mechanical torque (both calculated by initialize_machine)
    
    # Extract governor states
    states = zeros(Float64, 3)

    # Populate
    states[FV_IDX] = τ_m_init
    states[FF_IDX] = τ_m_init
    states[ET_IDX] = τ_m_init

    # Return things
    return states
end

# Update governor states
function update_gov_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    omega::Float64,
    gov::GasTG
)
    # Extract governor states
    fv = states[XG1_IDX]
    ff = states[XG2_IDX]
    et = states[XG3_IDX]

    # Calculate inverse droop (avoid division by zero) for numerical imstabilities faced when running the simulation
    inv_R = gov.R < eps() ? 0.0 : 1.0 / gov.R

    # Printing this value to make sure it is feasile
    println(inv_R)

    # Speed governor input calculation
    # TODO: replace 1.0 with \omega ref if needed.
    speed_error = omega - 1.0
    reference_power = gov.P_ref
    speed_governing = inv_R * speed_error

    # Temperature control limiter
    temperature_limit = gov.AT + gov.KT * (gov.AT - xg3)

    # Take minimum of power reference and temperature limit (LV logic gate)
    governor_input = min(reference_power - speed_governing, temperature_limit)

    # Process through low pass filters with anti-windup
    ff_sat, dff_dt = low_pass_nonwindup(governor_input, ff, 1.0, gov.T1, gov.V_min, gov.V_max)

    dfv_dt = low_pass(ff_sat, fv, 1.0, gov.T2)
    det_dt = low_pass(fv, et, 1.0, gov.T3)

    # Update derivatives
    derivatives[XG1_IDX] = dff_dt
    derivatives[XG2_IDX] = dfv_dt
    derivatives[XG3_IDX] = det_dt

    # Calculate mechanical power output
    # TODO: Check this since D_turb is 0.0
    Pm = fv - gov.D_turb * speed_error

    # Calculate mechanical torque (P = ω * τ, so τ = P/ω)
    τm = Pm / omega

    return τm
end

end # module