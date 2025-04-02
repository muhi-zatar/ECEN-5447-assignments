# Defining the module
module EXST1Model

# Exporting the needed functions
export EXST1, initialize_avr, update_avr_states!

# AVR state indices
const EFD_IDX = 1   # Field voltage (output of amplifier)
const VS_IDX = 2    # Sensed terminal voltage
const VLL_IDX = 3   # Lead-lag output
const VF_IDX = 4    # Feedback signal

# IEEE Type ST1 excitation system model structure
mutable struct EXST1
    # AVR parameters
    TR::Float64      # Voltage measurement time constant (s)
    TB::Float64      # Lead-lag denominator time constant (s)
    TC::Float64      # Lead-lag numerator time constant (s)
    KF::Float64      # Feedback gain (pu)
    TF::Float64      # Feedback time constant (s)
    KA::Float64      # Amplifier gain (pu)
    TA::Float64      # Amplifier time constant (s)
    KC::Float64      # Rectifier loading factor related to commutating reactance (pu)
    V_ref::Float64   # Reference Voltage

    # Constructor with default values
    function EXST1(;
        TR=0.01,
        TB=20.0,
        TC=10.0,
        KF=0.0,
        TF=0.1,
        KA=200.0,
        TA=0.1,
        KC=0.0,
        V_ref=1.0
    )
        return new(TR, TB, TC, KF, TF, KA, TA, KC, V_ref)
    end
end

# Initialize AVR states for proper steady-state equilibrium
function initialize_avr(avr::EXST1, Vt_init::Float64, Efd_init::Float64)
    # Initialize state vector
    states = zeros(Float64, 4)

    # Terminal voltage sensing
    states[VS_IDX] = Vt_init

    # For defd_dt = 0: KA*vll = efd
    states[VLL_IDX] = Efd_init / avr.KA

    # For dvll_dt = 0: (verr * (1.0 + TC/TB) - vll) = 0
    # Solving for verr: verr = vll / (1.0 + TC/TB)
    steady_state_verr = states[VLL_IDX] / (1.0 + avr.TC / avr.TB)

    # For dvf_dt = 0: vf = 0 (if defd_dt = 0)
    states[VF_IDX] = 0.0

    # Set Vref to achieve the required steady-state error
    # verr = Vref - vs - vf => Vref = verr + vs + vf
    avr.V_ref = steady_state_verr + states[VS_IDX] + states[VF_IDX]
    println("V_ref value is: $(avr.V_ref)")
    # Field voltage
    states[EFD_IDX] = Efd_init

    return states
end

# Update AVR states using state space formulation with stability focus
function update_avr_states!(
    states::AbstractVector{Float64},
    derivatives::AbstractVector{Float64},
    Vt::Float64,
    avr::EXST1,
)
    # Extract state variables
    efd = states[EFD_IDX]
    vs = states[VS_IDX]
    vll = states[VLL_IDX]
    vf = states[VF_IDX]

    # Calculate voltage error
    # VERR = Vref - VS - VF
    verr = avr.V_ref - vs - vf

    # State equations

    # Measured terminal voltage 
    # With time constant TR
    dvs_dt = (Vt - vs) / avr.TR

    # Lead-lag compensator 
    # Transfer function (1+sTC)/(1+sTB)
    dvll_dt = ((verr * (1.0 + avr.TC / avr.TB)) - vll) / avr.TB

    # Amplifier 
    # With gain KA and time constant TA
    defd_dt = (avr.KA * vll - efd) / avr.TA

    # Feedback path 
    # Calculate from the block diagram
    dvf_dt = (avr.KF * (avr.KA * vll - efd) / avr.TA - vf) / avr.TF

    # Update derivatives vector
    derivatives[EFD_IDX] = defd_dt
    derivatives[VS_IDX] = dvs_dt
    derivatives[VLL_IDX] = dvll_dt
    derivatives[VF_IDX] = dvf_dt

    # Output is the field voltage
    return efd
end

end # module