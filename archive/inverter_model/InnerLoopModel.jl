#==========================================================

 Implements the equations in PowerSimulationDyanmics.jl documentation
 This implements the equations:
 (1a) ξ̇_d = v_d_vi_ref - v_d
 (1b) ξ̇_q = v_q_vi_ref - v_q
 (1c) γ̇_d = i_d_cv_ref - i_d_cv
 (1d) γ̇_q = i_q_cv_ref - i_q_cv
 (1e) ϕ̇_d = ω_ad * (v_d - ϕ_d)
 (1f) ϕ̇_q = ω_ad * (v_q - ϕ_q)
 (1g) v_d_vi_ref = v_olc_ref - r_v * i_d + ω_olc * l_v * i_q
 (1h) v_q_vi_ref = -r_v * i_q - ω_olc * l_v * i_d
 (1i) i_d_cv_ref = k_pv * (v_d_vi_ref - v_d) + k_iv * ξ_d - c_f * ω_olc * v_q + k_ffi * i_d
 (1j) i_q_cv_ref = k_pv * (v_q_vi_ref - v_q) + k_iv * ξ_q + c_f * ω_olc * v_d + k_ffi * i_q
 (1k) v_d_ref-signal = k_pc * (i_d_cv_ref - i_d_cv) + k_ic * γ_d - ω_olc * l_f * i_q_cv + k_ffv * v_d - k_ad * (v_d - ϕ_d)
 (1l) v_q_ref-signal = k_pc * (i_q_cv_ref - i_q_cv) + k_ic * γ_q + ω_olc * l_f * i_d_cv + k_ffv * v_q - k_ad * (v_q - ϕ_q)
===========================================================#

module InnerLoopModel

# Exporting the needed functions
export InnerLoop, initialize_innerloop, update_innerloop_states!, XI_D_IDX, XI_Q_IDX, GAMMA_D_IDX,
    GAMMA_Q_IDX, PHI_D_IDX, PHI_Q_IDX

# Imports
using NLsolve

# Inner Loop state indices
const XI_D_IDX = 1
const XI_Q_IDX = 2
const GAMMA_D_IDX = 3
const GAMMA_Q_IDX = 4
const PHI_D_IDX = 5
const PHI_Q_IDX = 6

# Tolerance for NLsolve
const STRICT_NLSOLVE_F_TOLERANCE = 1e-8


# Inner Loop model structure
mutable struct InnerLoop
    # Current controller parameters
    k_pv::Float64     # Voltage controller proportional gain
    k_iv::Float64     # Voltage controller integral gain
    k_ffv::Float64    # Binary variable enabling the voltage feed-forward in output of current controllers
    r_v::Float64      # Virtual resistance in pu
    l_v::Float64      # Virtual inductance in pu
    k_pc::Float64     # Current controller proportional gain
    k_ic::Float64     # Current controller integral gain
    k_ffi::Float64    # Binary variable enabling current feed-forward in output of current controllers
    ω_ad::Float64     # Active damping low pass filter cut-off frequency
    k_ad::Float64     # Active damping gain
    # Constructor with default values
    function InnerLoop(;
        k_pv=0.59,
        k_iv=736.0,
        k_ffv=0.0,
        r_v=0.0,
        l_v=0.2,
        k_pc=1.27,
        k_ic=14.3,
        k_ffi=0.0,
        ω_ad=50.0,
        k_ad=0.2
    )
        return new(k_pv, k_iv, k_ffv, r_v, l_v, k_pc, k_ic, k_ffi, ω_ad, k_ad)
    end
end

# Initialize Inner Loop states
function initialize_innerloop(
    innerloop::InnerLoop,
    Id_inv::Float64,            # Inverter-side filter inductor current (in the filter/network SRF)
    Iq_inv::Float64,            # Inverter-side filter inductor current (in the filter/network SRF)
    Vd_flt::Float64,            # Filter capacitor voltage (in the filter/network SRF)
    Vq_flt::Float64,            # Filter capacitor voltage (in the filter/network SRF)
    Id_grd::Float64,            # Grid-side filter inductor current (in the filter/network SRF)
    Iq_grd::Float64,            # Grid-side filter inductor current (in the filter/network SRF)
    δθ_olc0::Float64,           # Outer-loop controller reference angle (initial guess)
    ω_olc0::Float64,            # Outer-loop controller frequency (initial guess)
    v_olc_ref0,                 # Outer-loop controller reference voltage (initial guess)
    V_dc::Float64,              # DC-link capacitor voltage
    Vd_inv::Float64,            # Inverter output voltage (from filter initialization)
    Vq_inv::Float64,            # Inverter output voltage (from filter initialization)
    c_f::Float64,               # Filter capacitance
    l_f::Float64                # Filter inductance
)

    # Use a non-linear solver function to find initial states
    function f!(out, x)
        δθ_olc = x[1]
        v_olc_ref = x[2]
        ξ_d = x[3]
        ξ_q = x[4]
        γ_d = x[5]
        γ_q = x[6]
        ϕ_d = x[7]
        ϕ_q = x[8]

        #----- Reference Frame Transformations -----#
        # Capacitor voltage from filter SRF to controller SRF
        V_dq = (Vd_flt + im * Vq_flt) * exp(-im * δθ_olc)
        v_d = real(V_dq)
        v_q = imag(V_dq)

        # Grid-side current from filter SRF to controller SRF
        I_dq = (Id_grd + im * Iq_grd) * exp(-im * δθ_olc)
        i_d = real(I_dq)
        i_q = imag(I_dq)

        # Inverter-side current from filter SRF to controller SRF
        I_dq_inv = (Id_inv + im * Iq_inv) * exp(-im * δθ_olc)
        i_d_cv = real(I_dq_inv)
        i_q_cv = imag(I_dq_inv)

        # Inverter output voltage from filter SRF to controller SRF
        V_dq_inv = (Vd_inv + im * Vq_inv) * exp(-im * δθ_olc)
        v_d_inv = real(V_dq_inv)
        v_q_inv = imag(V_dq_inv)

        #----- Algebraic Equations -----# 
        # Virtual impedance reference voltage calculations (1g and 1h)
        v_d_vi_ref = v_olc_ref - innerloop.r_v * i_d + ω_olc0 * innerloop.l_v * i_q
        v_q_vi_ref = -1 * innerloop.r_v * i_q - ω_olc0 * innerloop.l_v * i_d

        # Current controller reference current calculations (1i and 1j)
        i_d_cv_ref = innerloop.k_pv * (v_d_vi_ref - v_d) + innerloop.k_iv * ξ_d - c_f * ω_olc0 * v_q + innerloop.k_ffi * i_d
        i_q_cv_ref = innerloop.k_pv * (v_q_vi_ref - v_q) + innerloop.k_iv * ξ_q + c_f * ω_olc0 * v_d + innerloop.k_ffi * i_q

        # Reference signal voltage calculations (1k and 1l)     TODO: Find out what these are used for / what they are
        v_d_refsignal = innerloop.k_pc * (i_d_cv_ref - i_d_cv) + innerloop.k_ic * γ_d - ω_olc0 * l_f * i_q_cv + innerloop.k_ffv * v_d - innerloop.k_ad * (v_d - ϕ_d)
        v_q_refsignal = innerloop.k_pc * (i_q_cv_ref - i_q_cv) + innerloop.k_ic * γ_q + ω_olc0 * l_f * i_d_cv + innerloop.k_ffv * v_q - innerloop.k_ad * (v_q - ϕ_q)

        #----- Residuals -----#
        out[1] = v_d_vi_ref - v_d
        out[2] = v_q_vi_ref - v_q
        out[3] = i_d_cv_ref - i_d_cv
        out[4] = i_q_cv_ref - i_q_cv
        out[5] = v_d - ϕ_d
        out[6] = v_q - ϕ_q
        out[7] = v_d_refsignal - v_d_inv
        out[8] = v_q_refsignal - v_q_inv
    end

    # Initial guess and solution
    x0 = [δθ_olc0, v_olc_ref0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sol = NLsolve.nlsolve(f!, x0, ftol=STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Inner loop initialization failed")
    else
        #----- Extract solution -----#
        x0 = sol.zero
        states = zeros(Float64, 6)

        # Updated OLC variables
        # TODO: Figure out how/why this needs to be updated. Shouldn't the outer loop's initial guess be correct?
        δθ_olc = x0[1]
        v_olc_ref = x0[2]

        # States
        states[XI_D_IDX] = x0[3]
        states[XI_Q_IDX] = x0[4]
        states[GAMMA_D_IDX] = x0[5]
        states[GAMMA_Q_IDX] = x0[6]
        states[PHI_D_IDX] = x0[7]
        states[PHI_Q_IDX] = x0[8]

        # Converter modulation
        m0_dq = (Vd_inv + im * Vq_inv) * exp(-im * δθ_olc) / V_dc
        m0_d = real(m0_dq)
        m0_q = imag(m0_dq)

    end

    return states, δθ_olc, v_olc_ref, m0_d, m0_q
end

# Update Inner Loop states
function update_innerloop_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    Id_inv::Float64,            # Inverter-side filter inductor current (in the filter/network SRF)
    Iq_inv::Float64,            # Inverter-side filter inductor current (in the filter/network SRF)
    Vd_flt::Float64,            # Filter capacitor voltage (in the filter/network SRF)
    Vq_flt::Float64,            # Filter capacitor voltage (in the filter/network SRF)
    Id_grd::Float64,            # Grid-side filter inductor current (in the filter/network SRF)
    Iq_grd::Float64,            # Grid-side filter inductor current (in the filter/network SRF)
    δθ_olc::Float64,            # Outer-loop controller reference angle
    ω_olc::Float64,             # Outer-loop controller frequency
    v_olc_ref::Float64,         # Outer-loop controller reference voltage
    c_f::Float64,               # Filter capacitance
    l_f::Float64,               # Filter inductance
    innerloop::InnerLoop
)
    # Extract state variables
    ξ_d = states[XI_D_IDX]
    ξ_q = states[XI_Q_IDX]
    γ_d = states[GAMMA_D_IDX]
    γ_q = states[GAMMA_Q_IDX]
    ϕ_d = states[PHI_D_IDX]
    ϕ_q = states[PHI_Q_IDX]

    #----- Reference Frame Transformations -----#
    # Capacitor voltage from filter SRF to controller SRF
    V_dq = (Vd_flt + im * Vq_flt) * exp(-im * δθ_olc)
    v_d = real(V_dq)
    v_q = imag(V_dq)

    # Grid-side current from filter SRF to controller SRF
    I_dq = (Id_grd + im * Iq_grd) * exp(-im * δθ_olc)
    i_d = real(I_dq)
    i_q = imag(I_dq)

    # Inverter-side current from filter SRF to controller SRF
    I_dq_inv = (Id_inv + im * Iq_inv) * exp(-im * δθ_olc)
    i_d_cv = real(I_dq_inv)
    i_q_cv = imag(I_dq_inv)

    #----- Algebraic Equations -----#
    # Virtual impedance reference voltage calculations (1g and 1h)
    v_d_vi_ref = v_olc_ref - innerloop.r_v * i_d + ω_olc * innerloop.l_v * i_q
    v_q_vi_ref = -1 * innerloop.r_v * i_q - ω_olc * innerloop.l_v * i_d

    # Current controller reference current calculations (1i and 1j)
    i_d_cv_ref = innerloop.k_pv * (v_d_vi_ref - v_d) + innerloop.k_iv * ξ_d - c_f * ω_olc * v_q + innerloop.k_ffi * i_d
    i_q_cv_ref = innerloop.k_pv * (v_q_vi_ref - v_q) + innerloop.k_iv * ξ_q + c_f * ω_olc * v_d + innerloop.k_ffi * i_q

    # Reference signal voltage calculations (1k and 1l)
    v_d_refsignal = innerloop.k_pc * (i_d_cv_ref - i_d_cv) + innerloop.k_ic * γ_d - ω_olc * l_f * i_q_cv + innerloop.k_ffv * v_d - innerloop.k_ad * (v_d - ϕ_d)
    v_q_refsignal = innerloop.k_pc * (i_q_cv_ref - i_q_cv) + innerloop.k_ic * γ_q + ω_olc * l_f * i_d_cv + innerloop.k_ffv * v_q - innerloop.k_ad * (v_q - ϕ_q)

    #----- Differential Equations -----#
    derivatives[XI_D_IDX] = v_d_vi_ref - v_d                    # (1a)
    derivatives[XI_D_IDX] = v_q_vi_ref - v_q                    # (1b)
    derivatives[GAMMA_D_IDX] = i_d_cv_ref - i_d_cv              # (1c)
    derivatives[GAMMA_Q_IDX] = i_q_cv_ref - i_q_cv              # (1c)
    derivatives[PHI_D_IDX] = innerloop.ω_ad * (v_d - ϕ_d)       # (1e)
    derivatives[PHI_Q_IDX] = innerloop.ω_ad * (v_q - ϕ_q)       # (1f)

    return v_d_refsignal, v_q_refsignal
end

end # InnerLoopModel module