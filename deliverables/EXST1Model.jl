module EXST1Model

export EXST1, initialize_avr, update_avr_states!

# AVR state indices
const VM_IDX = 1
const VRLL_IDX = 2
const VR_IDX = 3
const VFB_IDX = 4

# EXST1 Automatic Voltage Regulator model structure
mutable struct EXST1
    # AVR parameters
    Tr::Float64      # Measurement time constant
    Vi_min::Float64  # Input limits minimum
    Vi_max::Float64  # Input limits maximum
    Tc::Float64      # Lead-lag numerator time constant
    Tb::Float64      # Lead-lag denominator time constant
    Ka::Float64      # Amplifier gain
    Ta::Float64      # Amplifier time constant
    Vr_min::Float64  # Regulator output minimum limit
    Vr_max::Float64  # Regulator output maximum limit
    Kc::Float64      # Rectifier loading factor
    Kf::Float64      # Stabilizer gain
    Tf::Float64      # Stabilizer time constant
    V_ref::Float64   # Reference voltage
    
    # Constructor with default values
    function EXST1(;
        Tr = 0.02,
        Vi_min = -999.0,
        Vi_max = 999.0,
        Tc = 1.0,
        Tb = 10.0,
        Ka = 200.0,
        Ta = 0.1,
        Vr_min = -5.0,
        Vr_max = 5.0,
        Kc = 0.0,
        Kf = 0.0,
        Tf = 1.0,
        V_ref = 1.0
    )
        return new(Tr, Vi_min, Vi_max, Tc, Tb, Ka, Ta, Vr_min, Vr_max, Kc, Kf, Tf, V_ref)
    end
end

# Helper functions for the AVR blocks
function clamp(value, min_val, max_val)
    return max(min(value, max_val), min_val)
end

function low_pass_mass_matrix(input, state, gain, time_constant)
    if time_constant > 0.0
        derivative = (gain * input - state) / time_constant
        return state, derivative
    else
        # If time constant is zero, output = gain * input (algebraic)
        return gain * input, 0.0
    end
end

function high_pass(input, state, gain, time_constant)
    if time_constant > 0.0
        output = gain * (input - state)
        derivative = input / time_constant - state / time_constant
        return output, derivative
    else
        # For zero time constant, output is zero
        return 0.0, 0.0
    end
end

function lead_lag_mass_matrix(input, state, gain, t_numerator, t_denominator)
    if t_denominator > 0.0
        output = (gain * t_numerator * input + state * (t_denominator - t_numerator)) / t_denominator
        derivative = (input - state) / t_denominator
        return output, derivative
    else
        # For zero denominator, output = gain * input
        return gain * input, 0.0
    end
end

# Initialize AVR states
function initialize_avr(avr::EXST1, V_terminal_magnitude, Xad_Ifd = 0.0)
    # For steady-state initialization, we assume all states are at equilibrium
    # Vm is the measured terminal voltage
    Vm = V_terminal_magnitude
    
    # At steady state, Vfb would be zero (output of high-pass filter)
    Vfb = 0.0
    
    # Calculate Vr at steady state - this would be the field voltage
    # adjusted by system constraints
    Vr = avr.V_ref
    
    # Vrll is the output of the lead-lag block
    Vrll = Vr / avr.Ka
    
    # Return initial states
    return [Vm, Vrll, Vr, Vfb]
end

# Update AVR states
function update_avr_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    V_terminal_magnitude::Float64,
    Vss::Float64,  # PSS output
    Ifd::Float64,  # Field current
    avr::EXST1
)
    # Extract AVR states
    Vm = states[VM_IDX]
    Vrll = states[VRLL_IDX]
    Vr = states[VR_IDX]
    Vfb = states[VFB_IDX]
    
    # Reference voltage with PSS correction
    V_ref = avr.V_ref + Vss
    
    # Calculate block outputs and derivatives
    _, dVm_dt = low_pass_mass_matrix(V_terminal_magnitude, Vm, 1.0, avr.Tr)
    
    # High pass filter for feedback
    y_hp, dVfb_dt = high_pass(Vr, Vfb, avr.Kf, avr.Tf)
    
    # Lead-lag compensator
    compensator_input = clamp(V_ref - Vm - y_hp, avr.Vi_min, avr.Vi_max)
    y_ll, dVrll_dt = lead_lag_mass_matrix(compensator_input, Vrll, 1.0, avr.Tc, avr.Tb)
    
    # Amplifier block
    _, dVr_dt = low_pass_mass_matrix(y_ll, Vr, avr.Ka, avr.Ta)
    
    # Update derivatives
    derivatives[VM_IDX] = dVm_dt
    derivatives[VRLL_IDX] = dVrll_dt
    derivatives[VR_IDX] = dVr_dt
    derivatives[VFB_IDX] = dVfb_dt
    
    # Calculate field voltage with rectifier loading adjustment
    Vf = clamp(Vr, V_terminal_magnitude * avr.Vr_min - avr.Kc * Ifd, 
                   V_terminal_magnitude * avr.Vr_max - avr.Kc * Ifd)
    
    return Vf
end

end # module