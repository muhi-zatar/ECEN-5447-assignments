module OuterLoopModel

# Exporting the needed functions
export OuterLoop, initialize_outerloop, update_outerloop_states!

# Outer Loop state indices
const P_CTRL_IDX = 1    # Active power control state
const Q_CTRL_IDX = 2    # Reactive power control state

# Outer Loop model structure
mutable struct OuterLoop
    # Active power droop parameters
    Rp::Float64       # Active power droop coefficient
    ωz::Float64       # Frequency setpoint (rad/s)
    
    # Reactive power control parameters
    Kpq::Float64      # Proportional gain for Q controller
    Kiq::Float64      # Integral gain for Q controller
    ωf::Float64       # Frequency for reactive power filter (rad/s)
    
    # Voltage controller parameters
    kpv::Float64      # Voltage controller proportional gain
    kiv::Float64      # Voltage controller integral gain
    kffv::Float64     # Binary enabling voltage feed-forward
    rv::Float64       # Virtual resistance (pu)
    lv::Float64       # Virtual inductance (pu)
    
    # Constructor with default values
    function OuterLoop(;
        Rp=0.05,
        ωz=2π*60,
        Kpq=2.0,
        Kiq=30.0,
        ωf=0.1322*π*50.0,
        kpv=0.59,
        kiv=736.0,
        kffv=0.0,
        rv=0.0,
        lv=0.2
    )
        return new(Rp, ωz, Kpq, Kiq, ωf, kpv, kiv, kffv, rv, lv)
    end
end

# Initialize Outer Loop states
function initialize_outerloop(outerloop::OuterLoop, p_init::Float64, q_init::Float64)
    # This function initializes the outer loop controller states
    
    return states
end

# Update Outer Loop states
function update_outerloop_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    p_meas::Float64,       # Measured active power
    q_meas::Float64,       # Measured reactive power
    vdc::Float64,          # DC link voltage
    pll_outputs::Vector{Float64},  # [δ_pll, ω_pll]
    p_ref::Float64,        # Active power reference
    q_ref::Float64,        # Reactive power reference
    outerloop::OuterLoop
)
    # Extract state variables
end

end # OuterLoopModel module