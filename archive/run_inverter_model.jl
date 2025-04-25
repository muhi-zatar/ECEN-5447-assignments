cd(@__DIR__)
using Pkg
Pkg.activate(".")

using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LinearAlgebra

include("NetworkModel.jl")
include("inverter_model/FilterModel.jl")
include("inverter_model/PLLModel.jl")
include("inverter_model/OuterLoopModel.jl")
include("inverter_model/InnerLoopModel.jl")
include("PowerFlowWrapper.jl")

using .NetworkModel
using .FilterModel
using .PLLModel
using .OuterLoopModel
using .InnerLoopModel
using .PowerFlowWrapper

# Definining some constants

# Network state indices
const BUS_INFINITE_BUS = 1
const BUS_MACHINE_MODEL = 2
const BUS_LOAD = 3
const NUM_STATES = 22
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
const I_1_D_IDX = 13
const I_1_Q_IDX = 14
const I_3_D_IDX = 15
const I_3_Q_IDX = 16
const I_B1_D_IDX = 17
const I_B1_Q_IDX = 18
const I_B2_D_IDX = 19
const I_B2_Q_IDX = 20
const I_B3_D_IDX = 21
const I_B3_Q_IDX = 22

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

# Parameters structure for the integrated model
mutable struct InverterParams
    network::ThreeBusNetwork
    filter::Filter
    pll::PLL
    outerloop::OuterLoop
    innerloop::InnerLoop
    ωsys::Float64
    Vdc::Float64
    network_idx::UnitRange{Int64}
    filter_idx::UnitRange{Int64}
    pll_idx::UnitRange{Int64}
    outerloop_idx::UnitRange{Int64}
    innerloop_idx::UnitRange{Int64}
end

function run_inverter_model(network_file)
    # Get power flow results
    V_sol, θ_sol, P_sol, Q_sol = get_powerflow_results(network_file)

    # Extract values for the machine bus
    V_mag = V_sol[BUS_MACHINE_MODEL]
    V_angle = θ_sol[BUS_MACHINE_MODEL]
    P = P_sol[BUS_MACHINE_MODEL]
    Q = Q_sol[BUS_MACHINE_MODEL]

    V_terminal = V_mag * exp(im * V_angle)
    println("Initial power flow results:")
    println("V_terminal = $V_mag p.u ∠ $V_angle rad")
    println("P = $P pu, Q = $Q pu")

    # Initialize components
    network = ThreeBusNetwork()
    filter = Filter()
    pll = PLL()
    outerloop = OuterLoop()
    innerloop = InnerLoop()

    # Initialize network states
    network_states, Id_grd, Iq_grd = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)

    # Extract terminal voltage at the inverter bus (in network DQ)
    Vd_grd = network_states[V_2_D_IDX]
    Vq_grd = network_states[V_2_Q_IDX]

    # System frequency (per unit)
    ωsys = 1.0

    # Initialize filter states
    v_term = [Vd_grd, Vq_grd]
    i_term = [Id_grd, Iq_grd]
    filter_states, Vd_inv, Vq_inv = initialize_filter(filter, v_term, i_term, ωsys)

    # Unpack initial filter states for other components to use
    Id_inv = filter_states[ID_INV]
    Iq_inv = filter_states[IQ_INV]
    Vd_flt = filter_states[VD_FLT]
    Vq_flt = filter_states[VQ_FLT]
    # NOTE: Id_grd, Iq_grd are filter states, but they've already been calculated by the network initialization

    # Convert filter capacitor voltage to RI for PLL initialization
    V_flt_ri = dq_ri(0.0) * [Vd_flt; Vq_flt]  # Convert DQ to RI

    # Convert grid-side current from the filter to RI for Outer Loop initialization
    I_flt_ri = dq_ri(0.0) * [Id_grd, Iq_grd]


    # Initialize PLL states with filter capacitor voltage
    pll_states = initialize_pll(pll, V_flt_ri)

    # Initialize Outer Loop states with lots of stuff
    outerloop_states, v_olc_ref0, δθ_olc0 = initialize_outerloop(outerloop, V_flt_ri, I_flt_ri)

    # Stand-in for DC-side model. TODO: Replace with actual model (?) V_dc is only used for getting m_dq for PWM model
    V_dc = 600.0

    # Initialize Inner Loop states with lots of stuff
    innerloop_states, δθ_olc, v_olc_ref, m0_d, m0_q = initialize_innerloop(innerloop, Id_inv, Iq_inv, Vd_flt, Vq_flt, Id_grd, Iq_grd, δθ_olc0, 1.0, v_olc_ref0, V_dc, Vd_inv, Vq_inv, filter.cf)

    # Combine all states
    states = vcat(network_states, filter_states, pll_states, outerloop_states, innerloop_states)

    # Define state indices
    network_idx = 1:22
    filter_idx = 23:28
    pll_idx = 29:31
    outerloop_idx = 32:34
    innerloop_idx = 35:40

    # Create parameter struct
    p = InverterParams(network, filter, pll, outerloop, innerloop, 1.0, V_dc, network_idx, filter_idx, pll_idx, outerloop_idx, innerloop_idx)

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

    mass_matrix = Diagonal(M_system)
    println("Mass matrix shape: $(size(mass_matrix))")

    # Print initial states summary
    println("Initial Line Currents:")
    println("I_12_D: $(network_states[I_12_D_IDX])")
    println("I_12_Q: $(network_states[I_12_Q_IDX])")

    println("Initial Filter States:")
    println("Id_inv: $(filter_states[ID_INV])")
    println("Vd_flt: $(filter_states[VD_FLT])")

    println("Initial PLL states:")
    println("vq_pll: $(pll_states[VQ_PLL_IDX])")
    println("epsilon_pll: $(pll_states[EPSILON_IDX])")
    println("theta_pll: $(pll_states[THETA_IDX])")

    println("Initial OuterLoop states:")
    println("δθ_olc0: $(outerloop_states[THETA_OLC])")
    println("P_M: $(outerloop_states[P_M])")
    println("P_M: $(outerloop_states[Q_M])")

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

    # System dynamics function
    function inverter_dynamics!(du, u, params, t)
        # Extract states for each component
        network_states = u[params.network_idx]
        filter_states = u[params.filter_idx]
        pll_states = u[params.pll_idx]
        outerloop_states = u[params.outerloop_idx]
        innerloop_states = u[params.innerloop_idx]

        # Make copies of states for Float64 compatibility
        network_states_f64 = convert.(Float64, network_states)
        filter_states_f64 = convert.(Float64, filter_states)
        pll_states_f64 = convert.(Float64, pll_states)
        outerloop_states_f64 = convert.(Float64, outerloop_states)
        innerloop_states_f64 = convert.(Float64, innerloop_states)

        # Arrays for derivatives
        du_network = similar(network_states_f64)
        du_filter = similar(filter_states_f64)
        du_pll = similar(pll_states_f64)
        du_outerloop = similar(outerloop_states_f64)
        du_innerloop = similar(innerloop_states_f64)

        # Extract terminal voltage and current from network at bus 2 (inverter bus)
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]
        v_grid = [v_2_d, v_2_q]

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
            params.ωsys,
            params.pll
        )

        # Extract PLL angle
        θ_pll = pll_states_f64[THETA_IDX]

        # 2. Update Outer Loop states
        δθ_olc, v_olc_ref, ω_olc = update_outerloop_states!(
            outerloop_states_f64,
            du_outerloop,
            v_flt_ri,
            i_grd_ri,
            params.ωsys,
            params.outerloop
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
            params.filter.cf,
            params.innerloop
        )

        # Prepare inverter voltage for filter update
        v_inv = [v_d_refsignal, v_q_refsignal]

        # 4. Update filter states
        update_filter_states!(
            filter_states_f64,
            du_filter,
            v_inv,
            v_grid,
            params.ωsys,
            params.filter
        )

        # Compute apparent power for the network
        P_terminal = v_2_d * i_grd_d + v_2_q * i_grd_q
        Q_terminal = v_2_q * i_grd_d - v_2_d * i_grd_q
        S_terminal = complex(P_terminal, Q_terminal)

        # 5. Update network states
        _, _, _, _, _ = update_network_states!(
            network_states_f64,
            du_network,
            S_terminal,
            params.network
        )

        # Copy derivatives to output
        du[params.network_idx] = du_network
        du[params.filter_idx] = du_filter
        du[params.pll_idx] = du_pll
        du[params.outerloop_idx] = du_outerloop
        du[params.innerloop_idx] = du_innerloop
        # print du values for debugging
        println("du_network: $du_network")
        println("du_filter: $du_filter")
        println("du_pll: $du_pll")
        println("du_outerloop: $du_outerloop")
        println("du_innerloop: $du_innerloop")
        # Print debugging info at integer time steps
        if abs(t - round(t)) < 0.00001
            println("t=$t: P=$(real(S_terminal)), Q=$(imag(S_terminal)), θ_pll=$(θ_pll), δθ_olc=$(δθ_olc)")
        end
    end

    inverter_system = ODEFunction(inverter_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 10.0)
    prob = ODEProblem(inverter_system, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [15.0]

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
        integrator.p.ωsys = 1.02  # 2% frequency increase

        # 3. Load change - uncommenting this will simulate a load change
        # integrator.p.network.Z_L *= 0.85  # 15% load increase
    end

    # Create a Callback function that represents the perturbation
    cb = DiscreteCallback(condition, affect!)

    # Run simulation
    sol = solve(prob, Rodas5(autodiff=false), dt=0.001, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    t = sol.t

    # Create derived quantities for plotting
    P_values = Float64[]
    Q_values = Float64[]
    v_inner_d_values = Float64[]
    v_inner_q_values = Float64[]

    for i in 1:length(t)
        # Get current state values
        network_states = sol[p.network_idx, i]
        filter_states = sol[p.filter_idx, i]
        innerloop_states = sol[p.innerloop_idx, i]

        # Terminal voltage and current
        v_2_d = network_states[V_2_D_IDX-p.network_idx[1]+1]
        v_2_q = network_states[V_2_Q_IDX-p.network_idx[1]+1]
        i_grd_d = filter_states[ID_GRD-p.filter_idx[1]+1]
        i_grd_q = filter_states[IQ_GRD-p.filter_idx[1]+1]

        # Calculate power
        P = v_2_d * i_grd_d + v_2_q * i_grd_q
        Q = v_2_q * i_grd_d - v_2_d * i_grd_q

        push!(P_values, P)
        push!(Q_values, Q)

        # Extract inner loop reference signals
        # Note: This is approximate as we don't store v_d_refsignal and v_q_refsignal directly
        # Using the states to reconstruct them would require re-running the inner loop update
        # Instead, we'll extract XI_D and XI_Q states as proxies
        push!(v_inner_d_values, innerloop_states[XI_D_IDX-p.innerloop_idx[1]+1])
        push!(v_inner_q_values, innerloop_states[XI_Q_IDX-p.innerloop_idx[1]+1])
    end

    # Create plots
    # Power
    p1 = plot(t, P_values,
        label="P", title="Active Power", linewidth=2)

    # Reactive Power
    p2 = plot(t, Q_values,
        label="Q", title="Reactive Power", linewidth=2)

    # Inner Control Voltage
    p3 = plot(t, v_inner_d_values,
        label="v_d", title="Inner Control Voltage", linewidth=2)
    plot!(p3, t, v_inner_q_values,
        label="v_q", linewidth=2)

    # Outer Control Angle
    p4 = plot(t, [sol[p.outerloop_idx[THETA_OLC], i] for i in 1:length(t)],
        label="δθ_olc", title="Outer Control Angle", linewidth=2)

    # PLL Angle
    p5 = plot(t, [sol[p.pll_idx[THETA_IDX], i] for i in 1:length(t)],
        label="θ_pll", title="PLL Angle", linewidth=2)

    # Network Voltages
    p6 = plot(t, [sol[p.network_idx[V_2_D_IDX], i] for i in 1:length(t)],
        label="Vd", title="Bus 2 Voltage", linewidth=2)
    plot!(p6, t, [sol[p.network_idx[V_2_Q_IDX], i] for i in 1:length(t)],
        label="Vq", linewidth=2)

    # Combine plots
    p_combined = plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(1200, 1800), left_margin=10mm)
    savefig(p_combined, "inverter_simulation_results.png")

    return sol
end

sol = run_inverter_model("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")