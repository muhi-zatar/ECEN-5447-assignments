# Defining the module
module SingleMassModel

# Exporting the required functions from this file
export SingleMass, initialize_shaft, update_shaft_states!

```
We are modeling the single mass model, so it only has 2 states.
```
# Shaft state indices
const DELTA_IDX = 1
const OMEGA_IDX = 2

# Single mass shaft model structure
mutable struct SingleMass
    # Shaft parameters
    H::Float64      # Inertia constant (s)
    D::Float64      # Damping coefficient (pu)
    system_base_frequency::Float64  # Base frequency (Hz)
    
    # Constructor with default values
    function SingleMass(;
        H = 5.0,
        D = 2.0,
        system_base_frequency = 60.0
    )
        return new(H, D, system_base_frequency)
    end
end

```
The module has two functions:
Initializing states, and updating states
This will help us in scalability
```

# Initialize shaft states
function initialize_shaft(shaft::SingleMass, initial_angle::Float64)
    # At steady state, omega is 1.0 pu
    omega = 1.0
    
    # Initial angle from power flow
    delta = initial_angle
    
    # Return initial states
    return [delta, omega]
end

# Update shaft states
function update_shaft_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    τe::Float64,     # Electromagnetic torque
    τm::Float64,     # Mechanical torque
    ω_sys::Float64,  # System frequency
    shaft::SingleMass
)
    # Extract shaft states
    delta = states[DELTA_IDX]
    omega = states[OMEGA_IDX]
    
    # Base frequency
    f0 = shaft.system_base_frequency
    
    ```
    Equations 15.5 in Milano's book
    ```

    # Compute derivativesy
    derivatives[DELTA_IDX] = 2.0 * π * f0 * (omega - ω_sys)
    
    # Speed derivative
    derivatives[OMEGA_IDX] = 1.0 / (2.0 * shaft.H) * (
        τm - τe - (shaft.D * (omega - 1.0)) / omega
    )
    
    return delta, omega
end

end # module