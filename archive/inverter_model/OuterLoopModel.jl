module OuterLoopModel

# Exporting the needed functions
export OuterLoop, initialize_outerloop, update_outerloop_states!

# Outer Loop state indices
const THETA_OLC = 1
const P_M = 2
const Q_M = 3

# Outer Loop model structure
mutable struct OuterLoop
    # Active power droop parameters
    Rp::Float64       # Active power droop coefficient
    ωz::Float64       # Frequency setpoint (rad/s)

    # Reactive power control parameters
    Kpq::Float64      # Proportional gain for Q controller
    Kiq::Float64      # Integral gain for Q controller
    ωf::Float64       # Frequency for reactive power filter (rad/s)

    # Reference values
    P_ref::Float64    # Active power reference
    Q_ref::Float64    # Reactive power reference
    v_ref::Float64    # Voltage reference
    ω_ref::Float64    # Frequency reference

    # Constructor with default values
    function OuterLoop(;
        Rp=0.05,
        ωz=2π * 60,
        Kpq=2.0,
        Kiq=30.0,
        ωf=0.1322 * π * 50.0,
        P_ref=0.0, # to be set in the initialization
        Q_ref=0.0, # to be set in the initialization
        ω_ref=1.0, # to be set in the initialization
        v_ref=1.0  # to be set in the initialization
    )
        return new(Rp, ωz, Kpq, Kiq, ωf, P_ref, Q_ref, v_ref, ω_ref)
    end
end

# Initialize Outer Loop states
function initialize_outerloop(outerloop::OuterLoop, V_filter_init::AbstractArray{Float64}, I_filter_init::AbstractArray{Float64})
    # This function initializes the outer loop controller states
    states = zeros(Float64, 3)

    vr_filter = V_filter_init[1]
    vi_filter = V_filter_init[2]
    ir_filter = I_filter_init[1]
    ii_filter = I_filter_init[2]


    # Calculating initial values for active and reactive
    p_e = vr_filter * ir_filter + vi_filter * ii_filter
    q_e = vi_filter * ir_filter - vr_filter * ii_filter

    # Inital value for θ_olc (unsure about this)
    θ_olc = atan(vi_filter / vr_filter)

    # Setting the reference values
    outerloop.P_ref = p_e
    outerloop.Q_ref = q_e
    outerloop.v_ref = sqrt(vr_filter^2 + vi_filter^2)

    # TODO: Check if these are correct:
    q_m = outerloop.Q_ref
    v_olc_ref = outerloop.v_ref

    # Settting the initial values
    states[THETA_OLC] = θ_olc
    states[P_M] = p_e
    states[Q_M] = q_e

    return states, v_olc_ref, θ_olc
end

# Update Outer Loop states
function update_outerloop_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    # vdc::Float64,
    V_filter::AbstractArray{Float64},
    I_filter::AbstractArray{Float64},
    ω_sys::Float64,
    outerloop::OuterLoop
)
    # Extracting the states
    θ_olc = states[THETA_OLC]
    p_m = states[P_M]
    q_m = states[Q_M]

    # Extracting Filter voltage and current
    vr_filter = V_filter[1]
    vi_filter = V_filter[2]
    ir_filter = I_filter[1]
    ii_filter = I_filter[2]

    # Calculating necessary parameters (TODO: extract 60.0 from base frequency (this is fine for now))
    Ωb = 2π * 60.0
    p_e = vr_filter * ir_filter + vi_filter * ii_filter
    q_e = vi_filter * ir_filter - vr_filter * ii_filter
    ω_olc = outerloop.ω_ref + outerloop.Rp * (outerloop.P_ref - p_m)

    # Not sure why this is needed, as it is not used in other equations
    v_olc_ref = outerloop.v_ref + outerloop.Kpq * (outerloop.Q_ref - q_m)

    # Calculating the derivatives as per the outer loop control equations in the documentation
    dθ_olc = Ωb * (ω_olc - ω_sys)
    dp_m = outerloop.ωz * (p_e - p_m)
    dq_m = outerloop.ωf * (q_e - q_m)

    # Updating the states
    derivatives[THETA_OLC] = dθ_olc
    derivatives[P_M] = dp_m
    derivatives[Q_M] = dq_m

    return θ_olc, v_olc_ref, ω_olc
end

end # OuterLoopModel module