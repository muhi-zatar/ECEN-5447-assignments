#==========================================================
 Reduced Order Phase-Locked Loop (PLL) Model
 Implemens the equations in PowerSimulationDyanmics.jl documentation
 This implements the equations:
 (2a) v̇q,pll = ωlp[vq,out - vq,pll]
 (2b) ε̇pll = vq,pll
 (2c) θ̇pll = Ωb·δωpll
 (2d) δωpll = 1.0 - ωsys + kp,pll·vq,pll + ki,pll·εpll
 (2e) ωpll = δωpll + ωsys
 (2f) vd,out + j·vq,out = (vr + j·vi)·e^(-j·θpll)
===========================================================#
module PLLModel
using NLsolve
# Exporting the needed functions
export PLL, initialize_pll, update_pll_states!

# PLL state indices
const VQ_PLL_IDX = 1    # PLL q-axis voltage state (vq,pll)
const EPSILON_IDX = 2   # PLL error integral state (εpll)
const THETA_IDX = 3     # PLL angle (θpll)

function ri_dq(δ::T) where {T<:Number}
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return Float64[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

# PLL model structure
mutable struct PLL
    # PLL parameters
    ωlp::Float64      # Low-pass filter cut-off frequency (rad/s)
    kppll::Float64    # PLL proportional gain
    kipll::Float64    # PLL integral gain
    Ωb::Float64 # Base frequency (rad/s)

    # Constructor with default values
    function PLL(;
        ωlp=1.322 * π * 50,
        kppll=2.0,
        kipll=20.0,
        Ωb=2π * 60
    )
        return new(ωlp, kppll, kipll, Ωb)
    end
end

# Initialize PLL states
function initialize_pll(pll::PLL, v_init::Vector{Float64})
    # This function initializes the PLL states
    # v_init: Initial voltage vector [vr, vi] in grid reference frame
    # theta_init: Initial phase angle (optional - if not provided, calculated from v_init)

    states = zeros(Float64, 3)

    # Extract grid voltage components (will be filter when intergrated)
    vr = v_init[1] # Real part of grid voltage
    vi = v_init[2] # Imaginary part of grid voltage

    # θ_init = atan(vi, vr) # Using atan2 equivalent for correct quadrant

    # println("Initial PLL angle (θ_init): $θ_init")
    # # # Calculate initial vq,out using equation (2f)
    # # # vd,out + j·vq,out = (vr + j·vi)·e^(-j·θpll)
    # # No need convert anything here, since from class we are initializing these by zeros.
    # # V_dq_pll = ri_dq(θ_init + pi / 2) * [vr; vi]
    # # vq_out = imag(V_dq_pll)
    θ_pll = atan(vi, vr)
    # println("Initial PLL angle (θ0_pll): $θ_pll")
    Vpll_q = 0.0
    ϵ_pll = 0.0

    states[VQ_PLL_IDX] = Vpll_q
    states[EPSILON_IDX] = ϵ_pll
    states[THETA_IDX] = θ_pll

    return states
end

# Update PLL states
function update_pll_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    v_out_ri::Vector{Float64},  # Converter output voltage [vr, vi] in grid/filter reference frame
    ωsys::Float64,            # System frequency in pu
    pll::PLL
)
    # Extract state variables
    vq_pll = states[VQ_PLL_IDX]     # vq,pll
    ϵ_pll = states[EPSILON_IDX] # εpll
    θ_pll = states[THETA_IDX]     # θpll

    # Extract converter output voltage components
    vr = v_out_ri[1] # Real part of converter output voltage
    vi = v_out_ri[2] # Imaginary part of converter output voltage

    # Calculate vd,out and vq,out using equation (2f)
    complex_v = (vr + im * vi) * exp(-im * θ_pll)
    vd_out = real(complex_v) # Not needed here
    vq_out = imag(complex_v)

    # Equation (2a): v̇q,pll = ωlp[vq,out - vq,pll]
    dvq_pll_dt = pll.ωlp * (vq_out - vq_pll)

    # Equation (2b): ε̇pll = vq,pll
    dϵ_pll_dt = vq_pll

    # Equation (2d): δωpll = 1.0 - ωsys + kp,pll·vq,pll + ki,pll·εpll
    δ_ω_pll = 1.0 - ωsys + pll.kppll * vq_pll + pll.kipll * ϵ_pll

    # Equation (2e): ωpll = δωpll + ωsys
    ω_pll = δ_ω_pll + ωsys

    # Equation (2c): θ̇pll = Ωb·δωpll
    dθ_pll_dt = pll.Ωb * δ_ω_pll

    # Update derivatives
    derivatives[VQ_PLL_IDX] = dvq_pll_dt
    derivatives[EPSILON_IDX] = dϵ_pll_dt
    derivatives[THETA_IDX] = dθ_pll_dt

    return
end

end # PLLModel module