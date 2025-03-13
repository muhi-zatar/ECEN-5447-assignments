# Defining the module
module GasTGModel

# Exporting the needed functions
export GasTG, initialize_gov, update_gov_states!

# Governor state indices
const XG1_IDX = 1
const XG2_IDX = 2
const XG3_IDX = 3

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
        R = 0.05,
        T1 = 0.4,
        T2 = 0.1,
        T3 = 3.0,
        D_turb = 0.0,
        AT = 1.0,
        KT = 2.0,
        V_min = 0.0,
        V_max = 1.0,
        P_ref = 1.0
    )
        return new(R, T1, T2, T3, D_turb, AT, KT, V_min, V_max, P_ref)
    end
end

```
The below equations are from the diagram in PowerWorld following link:
https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20GAST_PTI%20and%20GASTD.htm
```


# Helper functions for the governor blocks
# When checked against Jose's implementation in the julia package, they are not the same.
# Bunch of helper functions for the AVR blocks
function low_pass(input, state, gain, time_constant)
    if time_constant > 0.0
        output = state
        derivative = (gain * input - state) / time_constant
        return output, derivative
    else
        # For zero time constant, output follows input
        return gain * input, 0.0
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

```
The module has two functions:
Initializing states, and updating states
This will help us in scalability
```

# Initialize governor states
function initialize_gov(gov::GasTG, initial_power::Float64)
    # For steady-state initialization at a given power level
    
    # At steady state
    # xg3 = xg2 = initial_power (assuming no frequency deviation)
    xg3 = initial_power
    xg2 = initial_power
    
    # xg1 should be set to match the steady-state power
    xg1 = initial_power
    
    # Return initial states
    return [xg1, xg2, xg3]
end

# Update governor states
function update_gov_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    omega::Float64,
    gov::GasTG
)
    # Extract governor states
    xg1 = states[XG1_IDX]
    xg2 = states[XG2_IDX]
    xg3 = states[XG3_IDX]
    
    # Calculate inverse droop (avoid division by zero) for numerical imstabilities faced when running the simulation
    inv_R = gov.R < eps() ? 0.0 : 1.0 / gov.R
    
    # Speed governor input calculation
    speed_error = omega - 1.0
    reference_power = gov.P_ref
    speed_governing = inv_R * speed_error
    
    # Temperature control limiter
    temperature_limit = gov.AT + gov.KT * (gov.AT - xg3)
    
    # Take minimum of power reference and temperature limit
    governor_input = min(reference_power - speed_governing, temperature_limit)
    
    # Process through low pass filters with anti-windup
    xg1_sat, dxg1_dt = low_pass_nonwindup(governor_input, xg1, 1.0, gov.T1, gov.V_min, gov.V_max)
    _, dxg2_dt = low_pass(xg1_sat, xg2, 1.0, gov.T2)
    _, dxg3_dt = low_pass(xg2, xg3, 1.0, gov.T3)
    
    # Update derivatives
    derivatives[XG1_IDX] = dxg1_dt
    derivatives[XG2_IDX] = dxg2_dt
    derivatives[XG3_IDX] = dxg3_dt
    
    # Calculate mechanical power output
    # TODO: Check this since D_turb is 0.0
    Pm = xg2 - gov.D_turb * speed_error
    
    # Calculate mechanical torque (P = ω * τ, so τ = P/ω)
    τm = Pm / omega
    
    return τm
end

end # module