# defining the module
module EXST1Model

# Exporting the objects and functions
export EXST1, initialize_avr, update_avr_states!

# AVR state indices
const VF_IDX = 1
const VT_IDX = 2
const VLL_IDX = 3
const VFB_IDX = 4

# EXST1 Automatic Voltage Regulator model structure
mutable struct EXST1
    # AVR parameters
    # As defined in the PowerWorld document
    Tr::Float64      # Filter time constant
    Vi_min::Float64  # Input limits minimum
    Vi_max::Float64  # Input limits maximum
    Tc::Float64      # Lag numerator time constant
    Tb::Float64      # Lag denominator time constant
    Ka::Float64      # Amplifier gain
    Ta::Float64      # Voltage Regulator time constant
    Vr_min::Float64  # Maximum Control element output
    Vr_max::Float64  # Minimum Control element output
    Kc::Float64      # Rectifier loading factor
    Kf::Float64      # Rate feedback gain
    Tf::Float64      # Rate feedback time constant
    V_ref::Float64   # Reference voltage
    V_ss::Float64     # Stabilizer Output

    # Constructor with given values in the project description
    function EXST1(;
        Tr=0.01,
        Vi_min=-5.0,
        Vi_max=5.0,
        Tc=10.0,
        Tb=20.0,
        Ka=200.0,
        Ta=0.1,
        Vr_min=0,
        Vr_max=6,
        Kc=0.0,
        Kf=0.0,
        Tf=0.1,
        V_ref=1.0,
        V_ss=0.0
    )
        return new(Tr, Vi_min, Vi_max, Tc, Tb, Ka, Ta, Vr_min, Vr_max, Kc, Kf, Tf, V_ref)
    end
end

function clamp(value, min_val, max_val)
    return max(min(value, max_val), min_val)
end

function low_pass_mass_matrix(input, state, gain, time_constant)
    if time_constant > 0.0
        derivative = (gain * input - state) / time_constant
        return gain * input, derivative
    else
        return gain * input, 0.0
    end
end



function high_pass(input, state, gain, time_constant)
    if time_constant > 0.0
        derivative = (input - state) / time_constant
        new_state = state + derivative
        output = gain * new_state
        return output, derivative
    else
        return 0.0, 0.0
    end
end

function lead_lag_mass_matrix(input, state, gain, t_numerator, t_denominator)
    if t_denominator > 0.0
        output = (gain * t_numerator * input + state * (t_denominator - t_numerator)) / t_denominator
        derivative = (input - state) / t_denominator
        return output, derivative
    else
        return gain * input, 0.0
    end
end



# END OF HELPER FUNCTIONS


# The module has two functions:
# Initializing states, and updating states
# This will help us in scalability


# This function will initialize the AVR states based on the steady-state terminal voltage and the
# steady-state field voltage (as calculated by initialize_machine)

# Initialize AVR states
function initialize_avr(avr::EXST1, V_terminal_magnitude, V_f_init)

    # Set V_ref to field voltage calculated from initial power flow solution (this was set to V_f_init, I think it is wrong)
    avr.V_ref = V_terminal_magnitude

    # Set V_ss
    avr.V_ss = V_f_init / avr.Ka

    # Define empty states vector to populate
    states = zeros(Float64, 4)

    # Extract AVR states
    states[VF_IDX] = V_f_init
    states[VT_IDX] = V_terminal_magnitude
    states[VLL_IDX] = V_f_init / avr.Ka
    states[VFB_IDX] = 0

    return states
end

function update_avr_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    V_terminal_magnitude::Float64,
    avr::EXST1
)

    # Steps for calculation:
    # Step 1:Sense the Terminal Voltage and Apply a Low-Pass Filter (Block 2 in the Diagram)
    # Step 2: Compute the Voltage Error # Summation Block 
    # Step 3: Apply Lead-Lag Compensation (Block 3 in the Diagram (Lead-Lag Filter))
    # Step 4: Amplifier Response ( Block 1 in the Diagram )
    # Step 5: Compute Field Voltage
    # step 6: Apply Rate Feedback for Damping
    # Step 7: update derivates
    # Extract AVR states
    Vf = states[VF_IDX]
    Vt = states[VT_IDX]
    Vll = states[VLL_IDX]
    Vfb = states[VFB_IDX]

    Vt_filtered, dVt_dt = low_pass_mass_matrix(V_terminal_magnitude, Vt, 1.0, avr.Tr)

    y_hp, dVfb_dt = high_pass(Vf, Vfb, avr.Kf, avr.Tf)

    # Added stabilizer
    V_err = avr.V_ref - Vt_filtered
    compensator_input = V_err + avr.V_ss - y_hp

    y_ll, dVll_dt = lead_lag_mass_matrix(compensator_input, Vll, 1.0, avr.Tc, avr.Tb)

    # From regulator
    _, dVf_dt = low_pass_mass_matrix(y_ll, Vf, avr.Ka, avr.Ta)

    # Without clamping
    # Vf_new = V_reg
    # dVf_dt = (Vf_new - Vf) / avr.Ta

    derivatives[VF_IDX] = dVf_dt
    derivatives[VT_IDX] = dVt_dt
    derivatives[VLL_IDX] = dVll_dt
    derivatives[VFB_IDX] = dVfb_dt

    return
end



end # module