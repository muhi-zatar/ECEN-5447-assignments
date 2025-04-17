module InnerLoopModel

# Exporting the needed functions
export InnerLoop, initialize_innerloop, update_innerloop_states!

# Inner Loop state indices
const ID_CTRL_IDX = 1    # d-axis current control state
const IQ_CTRL_IDX = 2    # q-axis current control state

# Inner Loop model structure
mutable struct InnerLoop
    # Current controller parameters
    kpc::Float64      # Current controller proportional gain
    kic::Float64      # Current controller integral gain
    kffi::Float64     # Binary enabling current feed-forward
    ωad::Float64      # Active damping low pass filter cut-off frequency
    kad::Float64      # Active damping gain
    
    # Constructor with default values
    function InnerLoop(;
        kpc=1.27,
        kic=14.3,
        kffi=0.0,
        ωad=50.0,
        kad=0.2
    )
        return new(kpc, kic, kffi, ωad, kad)
    end
end

# Initialize Inner Loop states
function initialize_innerloop(innerloop::InnerLoop, id_init::Float64, iq_init::Float64)
    # This function initializes the inner loop controller states
    
    return states
end

# Update Inner Loop states
function update_innerloop_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    id_meas::Float64,      # Measured d-axis current
    iq_meas::Float64,      # Measured q-axis current
    id_ref::Float64,       # Reference d-axis current
    iq_ref::Float64,       # Reference q-axis current
    outerloop_outputs::Vector{Float64},  # [vdc, δ_olc, q_control]
    v_meas::Vector{Float64},  # Measured voltage [vd, vq]
    innerloop::InnerLoop
)
    # Extract state variables

end

end # InnerLoopModel module