# Ensure script runs in the correct environment
cd(@__DIR__)
using Pkg
Pkg.activate(".")

using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LinearAlgebra

include("NetworkModel_Multimachine.jl")
include("inverter_model/FilterModel.jl")
include("inverter_model/PLLModel.jl")
include("inverter_model/OuterLoopModel.jl")
include("inverter_model/InnerLoopModel.jl")
include("EXST1Model.jl")
include("GasTGModel.jl")
include("SauerPaiMachineModel.jl")
include("PowerFlowWrapper.jl")

using .NetworkModel_Multimachine
using .FilterModel
using .PLLModel
using .OuterLoopModel
using .InnerLoopModel
using .EXST1Model
using .GasTGModel
using .SauerPaiMachineModel
using .PowerFlowWrapper

# Definining some constants
const BUS_MACHINE_MODEL = 1
const BUS_CONVERTER_MODEL = 2
const BUS_LOAD = 3
const NUM_STATES = 20 # Updated to 20 as requested
const I_12_D_IDX = 1
const I_12_Q_IDX = 2
const I_13_D_IDX = 3
const I_13_Q_IDX = 4
const I_23_D_IDX = 5
const I_23_Q_IDX = 6
const V_1_D_IDX = 7
const V_1_Q_IDX = 8
const V_2_D_IDX = 9
const V_2_Q_IDX = 10
const V_3_D_IDX = 11
const V_3_Q_IDX = 12
const I_3_D_IDX = 13
const I_3_Q_IDX = 14
const I_B1_D_IDX = 15
const I_B1_Q_IDX = 16
const I_B2_D_IDX = 17
const I_B2_Q_IDX = 18
const I_B3_D_IDX = 19
const I_B3_Q_IDX = 20

# Filter state indices
const ID_INV = 1       # Inverter-side inductor current (real component)
const IQ_INV = 2       # Inverter-side inductor current (imag component)
const VD_FLT = 3       # Filter capacitor voltage (real component)
const VQ_FLT = 4       # Filter capacitor voltage (imag component)
const ID_GRD = 5       # Grid-side inductor current (real component)
const IQ_GRD = 6       # Grid-side inductor current (imag component)

# PLL state indices
const VQ_PLL_IDX = 1   # PLL q-axis voltage state (vq,pll)
const EPSILON_IDX = 2  # PLL error integral state (εpll)
const THETA_IDX = 3    # PLL angle (θpll)

# Outer Loop state indices
const THETA_OLC = 1
const P_M = 2
const Q_M = 3

# Inner loop state indices are exported by the module

# Machine state indices
const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6
const PSI_D = 7
const PSI_Q = 8

# AVR state indices
const EFD_IDX = 1   # Field voltage (output of amplifier)
const VS_IDX = 2    # Sensed terminal voltage
const VLL_IDX = 3   # Lead-lag output
const VF_IDX = 4    # Feedback signal

# Governor state indices
const FV_IDX = 1                # Fuel value
const FF_IDX = 2                # Fuel flow
const ET_IDX = 3                # Exhaust temp

# Matrix transformation functions
function ri_dq(δ::T) where {T<:Number}
    return Float64[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

function dq_ri(δ::T) where {T<:Number}
    return Float64[
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

# Function to transform from one dq reference frame to another
function dq_to_dq(v_dq::Vector{Float64}, δ_from::Float64, δ_to::Float64)
    # First convert to ri
    v_ri = dq_ri(δ_from) * v_dq
    # Then convert to the new dq frame
    return ri_dq(δ_to) * v_ri
end

# Parameters structure for the integrated model
mutable struct MultiMachineParams
    network::ThreeBusNetwork
    filter::Filter
    pll::PLL
    outerloop::OuterLoop
    innerloop::InnerLoop
    Vdc::Float64
    machine::SauerPaiMachine
    avr::EXST1
    governor::GasTG
    network_idx::UnitRange{Int64}
    filter_idx::UnitRange{Int64}
    pll_idx::UnitRange{Int64}
    outerloop_idx::UnitRange{Int64}
    innerloop_idx::UnitRange{Int64}
    machine_idx::UnitRange{Int64}
    avr_idx::UnitRange{Int64}
    gov_idx::UnitRange{Int64}
    ωsys::Float64  # Adding system frequency as a parameter
end

function run_multimachine_model(network_file)
    # Get power flow results
    V_sol, θ_sol, P_sol, Q_sol = get_powerflow_results(network_file)

    # Extract values for the machine bus
    V_mag_machine = V_sol[BUS_MACHINE_MODEL]
    V_angle_machine = θ_sol[BUS_MACHINE_MODEL]
    P_machine = P_sol[BUS_MACHINE_MODEL]
    Q_machine = Q_sol[BUS_MACHINE_MODEL]

    V_terminal_machine_init = V_mag_machine * exp(im * V_angle_machine)
    println("Initial power flow results for Machine bus:")
    println("V_terminal_machine = $V_mag_machine p.u ∠ $V_angle_machine rad")
    println("P = $P_machine pu, Q = $Q_machine pu")

    # Extract values for the converter bus
    V_mag_converter = V_sol[BUS_CONVERTER_MODEL]
    V_angle_converter = θ_sol[BUS_CONVERTER_MODEL]
    P_converter = P_sol[BUS_CONVERTER_MODEL]
    Q_converter = Q_sol[BUS_CONVERTER_MODEL]

    V_terminal_converter_init = V_mag_converter * exp(im * V_angle_converter)
    println("Initial power flow results for Converter bus:")
    println("V_terminal = $V_mag_converter p.u ∠ $V_angle_converter rad")
    println("P = $P_converter pu, Q = $Q_converter pu")

    # Initialize components
    network = ThreeBusNetwork()
    filter = Filter()
    pll = PLL()
    outerloop = OuterLoop()
    innerloop = InnerLoop()
    machine = SauerPaiMachine(
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

    avr = EXST1(
        V_ref=V_mag_machine  # Set reference voltage to initial voltage
    )

    governor = GasTG(
        P_ref=P_machine  # Set reference power to initial power
    )

    # System frequency (per unit)
    ωsys = 1.0

    # MODIFIED: Initialize machine first to get the reference angle
    machine_states, Vf_init, τ_m_init = initialize_machine(machine, 
                               V_terminal_machine_init, V_angle_machine, 
                               P_machine, Q_machine)
    
    # Extract initial machine angle for reference frame
    δ_machine_init = machine_states[DELTA]
    
    # Initialize network states with machine angle as reference
    network_states, Id_machine, Iq_machine, Id_grd, Iq_grd = initialize_network(
        network, V_sol, θ_sol, P_sol, Q_sol, δ_machine_init)

    # Debugging
    println("Initial Machine Angle (δ): $(δ_machine_init * 180/π) degrees")
    println("Initial Line Currents in machine reference frame:")
    println("I_12_D: $(network_states[I_12_D_IDX])")
    println("I_12_Q: $(network_states[I_12_Q_IDX])")
    println("I_13_D: $(network_states[I_13_D_IDX])")
    println("I_13_Q: $(network_states[I_13_Q_IDX])")
    println("I_23_D: $(network_states[I_23_D_IDX])")
    println("I_23_Q: $(network_states[I_23_Q_IDX])")

    # Extract terminal voltage at the inverter bus (in machine DQ)
    Vd_grd_machine = network_states[V_1_D_IDX]
    Vq_grd_machine = network_states[V_1_Q_IDX]

    # MODIFIED: Transform from machine reference to system reference for filter
    # First to RI coordinates
    v_term_ri = dq_ri(δ_machine_init) * [Vd_grd_machine; Vq_grd_machine]
    
    # Then to system reference frame (60Hz) for filter
    v_term_sys = ri_dq(0.0) * v_term_ri
    v_term = [v_term_sys[1], v_term_sys[2]]
    
    # Similarly for current
    i_term_ri = dq_ri(δ_machine_init) * [Id_grd; Iq_grd]
    i_term_sys = ri_dq(0.0) * i_term_ri
    i_term = [i_term_sys[1], i_term_sys[2]]

    # Initialize filter states with transformed values
    filter_states, Vd_inv, Vq_inv = initialize_filter(filter, v_term, i_term, ωsys)

    # Unpack initial filter states for other components to use
    Id_inv = filter_states[ID_INV]
    Iq_inv = filter_states[IQ_INV]
    Vd_flt = filter_states[VD_FLT]
    Vq_flt = filter_states[VQ_FLT]
    # NOTE: Id_grd, Iq_grd are filter states in system reference frame

    # Convert filter capacitor voltage to RI for PLL initialization
    V_flt_ri = dq_ri(0.0) * [Vd_flt; Vq_flt]

    # Convert grid-side current from the filter to RI for Outer Loop initialization
    I_flt_ri = dq_ri(0.0) * [filter_states[ID_GRD]; filter_states[IQ_GRD]]

    # Debugging
    println("Initial Filter States in system reference frame:")
    println("Id_inv: $(filter_states[ID_INV])")
    println("Iq_inv: $(filter_states[IQ_INV])")
    println("Vd_flt: $(filter_states[VD_FLT])")
    println("Vq_flt: $(filter_states[VQ_FLT])")
    println("Id_grd: $(filter_states[ID_GRD])")
    println("Iq_grd: $(filter_states[IQ_GRD])")
    println("Vr_flt: $(V_flt_ri[1])")
    println("Vi_flt: $(V_flt_ri[2])")
    println("Ir_flt: $(I_flt_ri[1])")
    println("Ii_flt: $(I_flt_ri[2])")

    println("Initial filter algebraic equations")
    println("Vd_inv: $(Vd_inv)")
    println("Vq_inv: $(Vq_inv)")

    # Initialize PLL states with filter capacitor voltage
    pll_states = initialize_pll(pll, V_flt_ri)

    # Debugging
    println("Initial PLL states:")
    println("vq_pll: $(pll_states[VQ_PLL_IDX])")
    println("epsilon_pll: $(pll_states[EPSILON_IDX])")
    println("theta_pll: $(pll_states[THETA_IDX])")

    # Initialize Outer Loop states with lots of stuff
    outerloop_states, v_olc_ref0, δθ_olc0 = initialize_outerloop(outerloop, V_flt_ri, I_flt_ri)

    # Debugging
    println("Initial OuterLoop states:")
    println("δθ_olc0: $(outerloop_states[THETA_OLC])")
    println("P_M: $(outerloop_states[P_M])")
    println("Q_M: $(outerloop_states[Q_M])")

    println("Initial OuterLoop algebraic variables:")
    println("v_olc_ref0: $(v_olc_ref0)")

    # Stand-in for DC-side model.
    V_dc = 1.0

    # Initialize Inner Loop states with lots of stuff
    innerloop_states, δθ_olc, v_olc_ref, m0_d, m0_q = initialize_innerloop(
        innerloop, Id_inv, Iq_inv, Vd_flt, Vq_flt, 
        filter_states[ID_GRD], filter_states[IQ_GRD], 
        δθ_olc0, 1.0, V_dc, Vd_inv, Vq_inv, filter.cf, filter.lf)

    println("Initial InnerLoop states:")
    println("ξ_d: $(innerloop_states[XI_D_IDX])")
    println("ξ_q: $(innerloop_states[XI_Q_IDX])")
    println("γ_d: $(innerloop_states[GAMMA_D_IDX])")
    println("γ_q: $(innerloop_states[GAMMA_Q_IDX])")
    println("ϕ_d: $(innerloop_states[PHI_D_IDX])")
    println("ϕ_q: $(innerloop_states[PHI_Q_IDX])")

    println("InnerLoop Initialization Has Updated:")
    println("δθ_olc: $(δθ_olc)")
    println("v_olc_ref: $(v_olc_ref)")

    println("InnerLoop Initialization Has Found:")
    println("m0_d: $(m0_d)")
    println("m0_q: $(m0_q)")

    # Override outer loop's initial guess at the angle and reference voltage
    outerloop_states[THETA_OLC] = δθ_olc
    set_V_ref(outerloop, v_olc_ref)

    # Initialize AVR and governor states
    avr_states = initialize_avr(avr, V_mag_machine, Vf_init)
    ω_init = 1.0
    governor_states = initialize_gov(governor, τ_m_init, ω_init)

    # Combine all states
    states = vcat(network_states, filter_states, pll_states, outerloop_states, innerloop_states, machine_states, avr_states, governor_states)
    
    # Print initial machine states
    println("Initial Machine States:")
    println("delta: $(machine_states[DELTA])")
    println("omega: $(machine_states[OMEGA])")
    println("Eq_p: $(machine_states[EQ_P])")
    println("Ed_p: $(machine_states[ED_P])")
    println("Psi_d_pp: $(machine_states[PSI_D_PP])")
    println("Psi_q_pp: $(machine_states[PSI_Q_PP])")
    println("Psi_d: $(machine_states[PSI_D])")
    println("Psi_q: $(machine_states[PSI_Q])")
    
    println("Initial AVR states:")
    println("Efd: $(avr_states[EFD_IDX])")
    println("Vs: $(avr_states[VS_IDX])")
    println("Vll: $(avr_states[VLL_IDX])")
    println("Vf: $(avr_states[VF_IDX])")
    
    println("Initial Governor states:")
    println("FV: $(governor_states[FV_IDX])")
    println("FF: $(governor_states[FF_IDX])")
    println("ET: $(governor_states[ET_IDX])")
    
    # Define state indices
    network_idx = 1:20  # Updated to 20 network states
    filter_idx = 21:26
    pll_idx = 27:29
    outerloop_idx = 30:32
    innerloop_idx = 33:38
    machine_idx = 39:46
    avr_idx = 47:50
    gov_idx = 51:53

    # Create parameter struct
    p = MultiMachineParams(
        network, filter, pll, outerloop, innerloop, V_dc, 
        machine, avr, governor, 
        network_idx, filter_idx, pll_idx, outerloop_idx, innerloop_idx, 
        machine_idx, avr_idx, gov_idx, 
        ωsys
    )

    # Build mass matrix
    M_system = zeros(Float64, length(states))

    if isa(network.M, Vector)
        M_system[network_idx] .= network.M
    else
        for i in 1:length(network_idx)
            M_system[network_idx[i]] = network.M[i, i]
        end
    end

    # Other states are differential
    M_system[filter_idx] .= 1.0
    M_system[pll_idx] .= 1.0
    M_system[outerloop_idx] .= 1.0
    M_system[innerloop_idx] .= 1.0
    M_system[machine_idx] .= 1.0
    M_system[avr_idx] .= 1.0
    M_system[gov_idx] .= 1.0

    mass_matrix = Diagonal(M_system)
    println("Mass matrix shape: $(size(mass_matrix))")

    # System dynamics function
    function multimachine_dynamics!(du, u, params, t)
        # ----- STEP 1: Rename variables for readability ----- # 
        # Extract states for each component
        network_states = u[params.network_idx]
        filter_states = u[params.filter_idx]
        pll_states = u[params.pll_idx]
        outerloop_states = u[params.outerloop_idx]
        innerloop_states = u[params.innerloop_idx]
        machine_states = u[params.machine_idx]
        avr_states = u[params.avr_idx]
        gov_states = u[params.gov_idx]

        # Extract parameters
        network = params.network
        filter = params.filter
        pll = params.pll
        outerloop = params.outerloop
        innerloop = params.innerloop
        machine = params.machine
        avr = params.avr
        governor = params.governor
        ωsys = params.ωsys

        # Make copies of states for Float64 compatibility
        network_states_f64 = convert.(Float64, network_states)
        filter_states_f64 = convert.(Float64, filter_states)
        pll_states_f64 = convert.(Float64, pll_states)
        outerloop_states_f64 = convert.(Float64, outerloop_states)
        innerloop_states_f64 = convert.(Float64, innerloop_states)
        machine_states_f64 = convert.(Float64, machine_states)
        avr_states_f64 = convert.(Float64, avr_states)
        gov_states_f64 = convert.(Float64, gov_states)

        # Arrays for derivatives
        du_network = similar(network_states_f64)
        du_filter = similar(filter_states_f64)
        du_pll = similar(pll_states_f64)
        du_outerloop = similar(outerloop_states_f64)
        du_innerloop = similar(innerloop_states_f64)
        du_machine = zeros(Float64, length(params.machine_idx))
        du_avr = zeros(Float64, length(params.avr_idx))
        du_gov = zeros(Float64, length(params.gov_idx))

        # ----- STEP 2: Update Machine States ----- # 
        # Extract machine state for reference frame transformation
        δ_machine = machine_states_f64[DELTA]
        ω_machine = machine_states_f64[OMEGA]
        
        # Grab terminal voltage from network at bus 2 (machine bus) in machine reference frame
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]
        
        # Convert to RI coordinates for machine update (note: already in machine reference frame)
        V_RI_machine = dq_ri(0.0) * [v_2_d; v_2_q]
        V_terminal_machine = complex(V_RI_machine[1], V_RI_machine[2])
        V_terminal_mag_machine = abs(V_terminal_machine)

        # Update AVR states
        update_avr_states!(avr_states_f64, du_avr, V_terminal_mag_machine, avr)
        
        # Grab field voltage from AVR
        efd = avr_states_f64[EFD_IDX]
        
        # Grab omega from machine
        ω = machine_states_f64[OMEGA]
        
        # Update Governor states
        τ_m = update_gov_states!(gov_states_f64, du_gov, ω, governor)
        
        # Update Machine states
        _, S_terminal_machine, _, _, _ = update_machine_states!(
            machine_states_f64, 
            du_machine, 
            V_terminal_machine, 
            efd, 
            τ_m, 
            machine
        )

        # ----- STEP 3: Update Inverter States ----- #
        # Extract terminal voltage from network at bus 1 (inverter bus) in machine reference frame
        v_1_d_machine = network_states_f64[V_1_D_IDX]
        v_1_q_machine = network_states_f64[V_1_Q_IDX]
        
        # Transform from machine reference to system reference for filter
        # First to RI coordinates
        v_grid_ri = dq_ri(δ_machine) * [v_1_d_machine; v_1_q_machine]
        
        # Then to system reference frame (60Hz) for filter
        v_grid_sys = ri_dq(0.0) * v_grid_ri
        v_grid = [v_grid_sys[1], v_grid_sys[2]]

        # Unpack filter states
        i_inv_d = filter_states_f64[ID_INV]
        i_inv_q = filter_states_f64[IQ_INV]
        v_flt_d = filter_states_f64[VD_FLT]
        v_flt_q = filter_states_f64[VQ_FLT]
        i_grd_d = filter_states_f64[ID_GRD]
        i_grd_q = filter_states_f64[IQ_GRD]

        # Convert filter capacitor voltage to rectangular coordinates for PLL
        v_flt_ri = dq_ri(0.0) * [v_flt_d; v_flt_q]

        # Convert grid-side filter current to rectangular coordinates for OuterLoop
        i_grd_ri = dq_ri(0.0) * [i_grd_d; i_grd_q]

        # 1. Update PLL states
        update_pll_states!(
            pll_states_f64,
            du_pll,
            v_flt_ri,
            ωsys,
            pll
        )

        # Extract PLL angle
        θ_pll = pll_states_f64[THETA_IDX]

        # 2. Update Outer Loop states
        δθ_olc, v_olc_ref, ω_olc = update_outerloop_states!(
            outerloop_states_f64,
            du_outerloop,
            v_flt_ri,
            i_grd_ri,
            ωsys,
            outerloop
        )

        # 3. Update Inner Loop states
        v_d_refsignal, v_q_refsignal = update_innerloop_states!(
            innerloop_states_f64,
            du_innerloop,
            i_inv_d,
            i_inv_q,
            v_flt_d,
            v_flt_q,
            i_grd_d,
            i_grd_q,
            δθ_olc,
            ω_olc,
            v_olc_ref,
            filter.cf,
            filter.lf,
            innerloop
        )

        # Prepare inverter voltage for filter update
        # Convert to network DQ to pass to network/filter
        v_inv_dq_network = (v_d_refsignal + im * v_q_refsignal) * exp(im * (δθ_olc + π / 2))

        # Pack into a vector for the filter
        v_inv_dq_network = [real(v_inv_dq_network), imag(v_inv_dq_network)]

        # 4. Update filter states
        update_filter_states!(
            filter_states_f64,
            du_filter,
            v_inv_dq_network,
            v_grid,
            ωsys,
            filter
        )

        # Compute apparent power for the network at inverter bus in system reference frame
        P_terminal_inverter_sys = v_grid[1] * i_grd_d + v_grid[2] * i_grd_q
        Q_terminal_inverter_sys = v_grid[2] * i_grd_d - v_grid[1] * i_grd_q
        S_terminal_inverter_sys = complex(P_terminal_inverter_sys, Q_terminal_inverter_sys)
        
        # Convert inverter current to machine reference frame for network update
        i_grd_ri = dq_ri(0.0) * [i_grd_d; i_grd_q]
        i_grd_machine = ri_dq(δ_machine) * i_grd_ri
        
        # Calculate power in machine reference frame
        P_terminal_inverter = v_1_d_machine * i_grd_machine[1] + v_1_q_machine * i_grd_machine[2]
        Q_terminal_inverter = v_1_q_machine * i_grd_machine[1] - v_1_d_machine * i_grd_machine[2]
        S_terminal_inverter = complex(P_terminal_inverter, Q_terminal_inverter)

        # ----- STEP 4: Update Network States ----- #
        # Now update network with both machine and inverter power injections
        # Pass machine angle and speed to the network update function
        _, _, _, _, _ = update_network_states!(
            network_states_f64,
            du_network,
            S_terminal_machine,  # Machine injection
            S_terminal_inverter, # Inverter injection in machine reference frame
            network,
            δ_machine,          # Pass machine angle
            ω_machine           # Pass machine speed
        )

        # Copy derivatives to output
        du[params.network_idx] = du_network
        du[params.filter_idx] = du_filter
        du[params.pll_idx] = du_pll
        du[params.outerloop_idx] = du_outerloop
        du[params.innerloop_idx] = du_innerloop
        du[params.machine_idx] = du_machine
        du[params.avr_idx] = du_avr
        du[params.gov_idx] = du_gov

        # println("du_network: $du_network")
        # println("du_filter: $du_filter")
        # println("du_pll: $du_pll")
        # println("du_outerloop: $du_outerloop")
        # println("du_innerloop: $du_innerloop")
        # println("du_machine: $du_machine")
        # println("du_avr: $du_avr")
        # println("du_gov: $du_gov")
        # exit()
        # Print debugging info at integer time steps
        if abs(t - round(t)) < 0.00001
            println("t=$t: Machine P=$(real(S_terminal_machine)), Q=$(imag(S_terminal_machine)), Inverter P=$(real(S_terminal_inverter)), Q=$(imag(S_terminal_inverter))")
            println("t=$t: Machine δ=$(δ_machine * 180/π)°, ω=$(ω_machine)")
        end
    end

    # Create ODEFunction
    multimachine_system = ODEFunction(multimachine_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 5.0)
    prob = ODEProblem(multimachine_system, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [1.0]

    # Define the condition for which to apply a perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end

    # Define the perturbation
    function affect!(integrator)
        # Choose one perturbation scenario

        # 1. Line Trip - uncommenting this will simulate a line trip
        # integrator.p.network.R_12 = 1e6
        # integrator.p.network.X_12 = 1e6
        # integrator.p.network.B_1 *= 0.5
        # integrator.p.network.B_2 *= 0.5

        # 2. Frequency change - uncommenting this will simulate a frequency change
        # integrator.p.ωsys = 1.02  # 2% frequency increase

        # 3. Load increase - uncommenting this will simulate a load increase
        # integrator.p.network.Z_L *= 0.85  # 15% load increase

        # 4. Load decrease - uncommenting this will simulate a load decrease
        integrator.p.network.Z_L *= 1.15
    end

    # Create a Callback function that represents the perturbation
    cb = DiscreteCallback(condition, affect!)

    # Run simulation
    sol = solve(prob, Rodas5(autodiff=false), dt=0.0001, adaptive=false, saveat=0.0001, callback=cb, tstops=perturb_times)

    t = sol.t

    # Allocate vectors for plotting
    line_12_current_values_r = Float64[]
    line_13_current_values_r = Float64[]
    line_23_current_values_r = Float64[]
    line_12_current_values_q = Float64[]
    line_13_current_values_q = Float64[]
# Allocate vectors for plotting
    line_23_current_values_q = Float64[]
    voltage_magnitude_values_1 = Float64[]
    voltage_magnitude_values_2 = Float64[]
    voltage_magnitude_values_3 = Float64[]
    bus_1_p = Float64[]
    bus_2_p = Float64[]
    bus_3_p = Float64[]
    bus_1_q = Float64[]
    bus_2_q = Float64[]
    bus_3_q = Float64[]
    ω_olc_values = Float64[]
    ω_pll_values = Float64[]
    machine_delta_values = Float64[]
    machine_omega_values = Float64[]
    avr_efd_values = Float64[]
    governor_torque_values = Float64[]

    # Add values to plotting vectors
    for i in 1:length(t)
        # Get current state values
        network_states = sol[p.network_idx, i]
        filter_states = sol[p.filter_idx, i]
        innerloop_states = sol[p.innerloop_idx, i]
        outerloop_states = sol[p.outerloop_idx, i]
        pll_states = sol[p.pll_idx, i]
        machine_states = sol[p.machine_idx, i]
        avr_states = sol[p.avr_idx, i]
        gov_states = sol[p.gov_idx, i]
        
        # Get machine angle for reference frame transformations
        δ_machine = machine_states[DELTA-p.machine_idx[1]+1]

        # Line current (already in machine reference frame)
        I_12_d = network_states[I_12_D_IDX-p.network_idx[1]+1]
        I_12_q = network_states[I_12_Q_IDX-p.network_idx[1]+1]
        I_13_d = network_states[I_13_D_IDX-p.network_idx[1]+1]
        I_13_q = network_states[I_13_Q_IDX-p.network_idx[1]+1]
        I_23_d = network_states[I_23_D_IDX-p.network_idx[1]+1]
        I_23_q = network_states[I_23_Q_IDX-p.network_idx[1]+1]

        # Terminal voltage in machine reference frame
        v_1_d = network_states[V_1_D_IDX-p.network_idx[1]+1]
        v_1_q = network_states[V_1_Q_IDX-p.network_idx[1]+1]
        v_3_d = network_states[V_3_D_IDX-p.network_idx[1]+1]
        v_3_q = network_states[V_3_Q_IDX-p.network_idx[1]+1]
        v_2_d = network_states[V_2_D_IDX-p.network_idx[1]+1]
        v_2_q = network_states[V_2_Q_IDX-p.network_idx[1]+1]
        
        # Filter grid-side current in system reference frame
        i_grd_d = filter_states[ID_GRD]
        i_grd_q = filter_states[IQ_GRD]
        
        # Convert filter current to machine reference frame for calculations
        i_grd_ri = dq_ri(0.0) * [i_grd_d; i_grd_q]
        i_grd_machine = ri_dq(δ_machine) * i_grd_ri
        
        # Machine states
        delta = machine_states[DELTA-p.machine_idx[1]+1]
        omega = machine_states[OMEGA-p.machine_idx[1]+1]
        
        # AVR and Governor
        efd = avr_states[EFD_IDX-p.avr_idx[1]+1]
        τ_m = governor.KT * (gov_states[FV_IDX-p.gov_idx[1]+1] - governor.Wf * gov_states[FF_IDX-p.gov_idx[1]+1])

        # Outer loop states
        p_m = outerloop_states[P_M-p.outerloop_idx[1]+1]
        ω_olc = outerloop.ω_ref + outerloop.Rp * (outerloop.P_ref - p_m)

        # PLL states
        vq_pll = pll_states[VQ_PLL_IDX-p.pll_idx[1]+1]
        ϵ_pll = pll_states[EPSILON_IDX-p.pll_idx[1]+1]
        δ_ω_pll = 1.0 - ωsys + pll.kppll * vq_pll + pll.kipll * ϵ_pll
        ω_pll = δ_ω_pll + ωsys

        # Convert line current from machine DQ to rectangular coordinates
        I_12_ri = dq_ri(0.0) * [I_12_d; I_12_q]
        I_13_ri = dq_ri(0.0) * [I_13_d; I_13_q]
        I_23_ri = dq_ri(0.0) * [I_23_d; I_23_q]

        # Calculate voltage magnitude
        voltage_magnitude_1 = sqrt(v_1_d^2 + v_1_q^2)
        voltage_magnitude_2 = sqrt(v_2_d^2 + v_2_q^2)
        voltage_magnitude_3 = sqrt(v_3_d^2 + v_3_q^2)

        # Calculate power at machine bus (already in machine reference frame)
        I_B2_d = network_states[I_B2_D_IDX-p.network_idx[1]+1]
        I_B2_q = network_states[I_B2_Q_IDX-p.network_idx[1]+1]
        P_2 = v_2_d * I_B2_d + v_2_q * I_B2_q
        Q_2 = v_2_q * I_B2_d - v_2_d * I_B2_q

        # Calculate power at inverter bus using machine reference frame
        P_1 = v_1_d * i_grd_machine[1] + v_1_q * i_grd_machine[2]
        Q_1 = v_1_q * i_grd_machine[1] - v_1_d * i_grd_machine[2]

        # Calculate power at load bus
        I_B3_d = network_states[I_B3_D_IDX-p.network_idx[1]+1]
        I_B3_q = network_states[I_B3_Q_IDX-p.network_idx[1]+1]
        P_3 = v_3_d * I_B3_d + v_3_q * I_B3_q
        Q_3 = v_3_q * I_B3_d - v_3_d * I_B3_q

        # Push values to plotting vectors
        push!(line_12_current_values_r, I_12_ri[1])
        push!(line_13_current_values_r, I_13_ri[1])
        push!(line_23_current_values_r, I_23_ri[1])
        push!(line_12_current_values_q, I_12_ri[2])
        push!(line_13_current_values_q, I_13_ri[2])
        push!(line_23_current_values_q, I_23_ri[2])
        push!(voltage_magnitude_values_1, voltage_magnitude_1)
        push!(voltage_magnitude_values_2, voltage_magnitude_2)
        push!(voltage_magnitude_values_3, voltage_magnitude_3)
        push!(bus_1_p, P_1)
        push!(bus_2_p, P_2)
        push!(bus_3_p, P_3)
        push!(bus_1_q, Q_1)
        push!(bus_2_q, Q_2)
        push!(bus_3_q, Q_3)
        push!(ω_olc_values, ω_olc)
        push!(ω_pll_values, ω_pll)
        push!(machine_delta_values, delta)
        push!(machine_omega_values, omega)
        push!(avr_efd_values, efd)
        push!(governor_torque_values, τ_m)
    end

    # Create plots
    # Active Power
    p1 = plot(t, bus_1_p,
    label="Inverter", title="Active Power", linewidth=2)
    plot!(p1, t, bus_2_p, label="Machine", linewidth=2)
    plot!(p1, t, bus_3_p, label="Load", linewidth=2)

    # Reactive Power
    p2 = plot(t, bus_1_q,
    label="Inverter", title="Reactive Power", linewidth=2)
    plot!(p2, t, bus_2_q, label="Machine", linewidth=2)
    plot!(p2, t, bus_3_q, label="Load", linewidth=2)

    # Outer Control Angle and PLL Angle in degrees
    p3 = plot(t, [sol[p.outerloop_idx[1]+THETA_OLC-1, i] * 180.0 / π for i in 1:length(t)],
    label="δθ_olc", title="Inverter Angles (Degrees)", linewidth=2)
    plot!(p3, t, [sol[p.pll_idx[1]+THETA_IDX-1, i] * 180.0 / π for i in 1:length(t)],
    label="θ_pll", linewidth=2)

    # Machine Rotor Angle and Speed
    p4 = plot(t, machine_delta_values * 180.0 / π,
    label="δ", title="Machine States", linewidth=2)
    plot!(p4, t, machine_omega_values, label="ω", linewidth=2, ylabel="Angle (deg) / Speed (pu)")

    # Plot Voltage Magnitude
    p5 = plot(t, voltage_magnitude_values_1,
    label="Inverter", title="Voltage Magnitude", linewidth=2)
    plot!(p5, t, voltage_magnitude_values_2, label="Machine", linewidth=2)
    plot!(p5, t, voltage_magnitude_values_3, label="Load", linewidth=2)

    # Plot Frequency
    p6 = plot(t, ω_pll_values .* 60.0, 
    label="PLL", linewidth=2, title="Frequency Comparison")
    plot!(p6, t, ω_olc_values .* 60.0, 
    label="Inv. Controller", linewidth=2)
    plot!(p6, t, machine_omega_values .* 60.0, 
    label="Machine", linewidth=2)

    # AVR and Governor
    p7 = plot(t, avr_efd_values,
    label="Field Voltage", linewidth=2, title="AVR Output")

    p8 = plot(t, governor_torque_values,
    label="Mech. Torque", linewidth=2, title="Governor Output")

    # Combine plots
    p_combined = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(4, 2), size=(1200, 1800), left_margin=10mm)
    savefig(p_combined, "multimachine_simulation_results.png")

    # Logging the states and writing to CSV
    open("multimachine_simulation_results.csv", "w") do io
        # Write the header
        write(io, "Time,P_1,Q_1,P_2,Q_2,P_3,Q_3,theta_pll,theta_olc,delta,omega,efd,torque,voltage_magnitude_1,voltage_magnitude_2,voltage_magnitude_3,I_12_r,I_12_i,I_13_r,I_13_i,I_23_r,I_23_i,omega_olc,omega_pll\n")
        # Write the data
        for i in 1:length(t)
            write(io, "$(t[i]),")
            write(io, "$(bus_1_p[i]),")
            write(io, "$(bus_1_q[i]),")
            write(io, "$(bus_2_p[i]),")
            write(io, "$(bus_2_q[i]),")
            write(io, "$(bus_3_p[i]),")
            write(io, "$(bus_3_q[i]),")
            write(io, "$(sol[p.pll_idx[1]+THETA_IDX-1, i] * 180.0 / π),")
            write(io, "$(sol[p.outerloop_idx[1]+THETA_OLC-1, i] * 180.0 / π),")
            write(io, "$(machine_delta_values[i] * 180.0 / π),")
            write(io, "$(machine_omega_values[i]),")
            write(io, "$(avr_efd_values[i]),")
            write(io, "$(governor_torque_values[i]),")
            write(io, "$(voltage_magnitude_values_1[i]),")
            write(io, "$(voltage_magnitude_values_2[i]),")
            write(io, "$(voltage_magnitude_values_3[i]),")
            write(io, "$(line_12_current_values_r[i]),")
            write(io, "$(line_12_current_values_q[i]),")
            write(io, "$(line_13_current_values_r[i]),")
            write(io, "$(line_13_current_values_q[i]),")
            write(io, "$(line_23_current_values_r[i]),")
            write(io, "$(line_23_current_values_q[i]),")
            write(io, "$(ω_olc_values[i]),")
            write(io, "$(ω_pll_values[i])\n")
        end
    end
    return sol
end

# Run the simulation
sol = run_multimachine_model("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")