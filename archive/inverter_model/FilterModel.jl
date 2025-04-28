module FilterModel

# Exporting the needed functions
export Filter, initialize_filter, update_filter_states!

# Filter state indices
const NUM_FLT_STATES = 6
const ID_INV = 1       # Inverter-side inductor current (real component)
const IQ_INV = 2       # Inverter-side inductor current (imag component)
const VD_FLT = 3       # Filter capacitor voltage (real component)
const VQ_FLT = 4       # Filter capacitor voltage (imag component)
const ID_GRD = 5       # Grid-side inductor current (real component)
const IQ_GRD = 6       # Grid-side inductor current (imag component)

# Filter model structure
"""
Notes on the filter
---------------------------------------------------------------------------------------------------
The filter is based on the model presented by D'Arco in the paper titled "A Virtual Synchronous 
Machine implementation for distributed control of power converters in SmartGrids."

In our model, we will implement the filter in the reference frame of the network, which will avoid
the need to transform the inverter voltage through two different reference frames.
"""
mutable struct Filter
    # Filter parameters
    lf::Float64         # Filter inductance (pu)
    rf::Float64         # Filter resistance (pu)
    cf::Float64         # Filter capacitance (pu)
    lg::Float64         # Grid inductance (pu)
    rg::Float64         # Grid resistance (pu)
    base_power::Float64 # Base power in MVA
    base_voltage::Float64 # Base voltage in kV
    Ωb::Float64 # Base frequency in rad/s

    # Constructor with default values
    function Filter(;
        lf=0.08,
        rf=0.003,
        cf=0.074,
        lg=0.2,
        rg=0.01,
        base_power=100.0,  # 100 MVA
        base_voltage=400.0, # Assuming 400 kV as base, adjust as needed
        Ωb=2π * 60   # 60 Hz system
    )
        return new(lf, rf, cf, lg, rg, base_power, base_voltage, Ωb)
    end
end

# Initialize filter states
"""
    Since the filter is in the same reference frame as the network, we will take advantage of the 
    fact that the network has already translated the power flow solution into d,q components for 
    the voltage and current at the terminal.
"""
function initialize_filter(flt::Filter, v_term::Vector{Float64}, i_term::Vector{Float64}, ωg::Float64)
    # Initialize filter states
    states = zeros(Float64, NUM_FLT_STATES)

    # Unpack the voltage components
    Vd_grd = v_term[1]
    Vq_grd = v_term[2]

    # Unpack the current components
    Id_grd = i_term[1]
    Iq_grd = i_term[2]

    # Start with the 5th and 6th differential equations
    Vd_flt = Vd_grd + flt.rg * Id_grd - ωg * flt.lg * Iq_grd
    Vq_flt = Vq_grd + flt.rg * Iq_grd + ωg * flt.lg * Id_grd

    # Then use the 3rd and 4th differential equations
    Id_inv = Id_grd - ωg * flt.cf * Vq_flt
    Iq_inv = Iq_grd + ωg * flt.cf * Vd_flt

    # Finally use the 1st and 2nd differential equations
    Vd_inv = Vd_flt + flt.rf * Id_inv - ωg * flt.lf * Iq_inv
    Vq_inv = Vq_flt + flt.rf * Iq_inv + ωg * flt.lf * Id_inv

    # Pack the states vector
    states[ID_INV] = Id_inv
    states[IQ_INV] = Iq_inv
    states[VD_FLT] = Vd_flt
    states[VQ_FLT] = Vq_flt
    states[ID_GRD] = Id_grd
    states[IQ_GRD] = Iq_grd

    # Return the states as well as the inverter voltage to serve as control points for the inner loop
    return states, Vd_inv, Vq_inv
end

# Update filter states using state space formulation
function update_filter_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    v_inv::Vector{Float64},     # Inverter voltage [vd, vq]
    v_grid::Vector{Float64},    # Grid voltage [vd, vq]
    ωg::Float64,             # Grid frequency
    flt::Filter
)
    # Extract state variables
    Id_inv = states[ID_INV]
    Iq_inv = states[IQ_INV]
    Vd_flt = states[VD_FLT]
    Vq_flt = states[VQ_FLT]
    Id_grd = states[ID_GRD]
    Iq_grd = states[IQ_GRD]

    # Unpack voltage vectors
    Vd_inv, Vq_inv = v_inv
    Vd_grd, Vi_grd = v_grid

    println("In update:")
    println("Vd_inv = $Vd_inv")
    println("Vq_inv = $Vq_inv")
    println("Vd_grd = $Vd_grd")
    println("Vq_grd = $Vq_grd")

    # State equations for the filter
    # Inverter-side inductor current
    derivatives[ID_INV] = (flt.Ωb / flt.lf) * (Vd_inv - Vd_flt - flt.rf * Id_inv + ωg * flt.lf * Iq_inv)
    derivatives[IQ_INV] = (flt.Ωb / flt.lf) * (Vq_inv - Vq_flt - flt.rf * Iq_inv - ωg * flt.lf * Id_inv)

    # Filter capacitor voltage
    derivatives[VD_FLT] = (flt.Ωb / flt.cf) * (Id_inv - Id_grd + ωg * flt.cf * Vq_flt)
    derivatives[VQ_FLT] = (flt.Ωb / flt.cf) * (Iq_inv - Iq_grd - ωg * flt.cf * Vd_flt)

    # Grid-side inductor current
    derivatives[ID_GRD] = (flt.Ωb / flt.lg) * (Vd_flt - Vd_grd - flt.rg * Id_grd + ωg * flt.lg * Iq_grd)
    derivatives[IQ_GRD] = (flt.Ωb / flt.lg) * (Vq_flt - Vi_grd - flt.rg * Iq_grd - ωg * flt.lg * Id_grd)
end

end # FilterModel module