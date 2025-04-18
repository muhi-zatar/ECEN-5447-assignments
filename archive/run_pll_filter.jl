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
include("PowerFlowWrapper.jl")

using .NetworkModel
using .FilterModel
using .PLLModel
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
mutable struct PLLFilterParams
    network::ThreeBusNetwork
    filter::Filter
    pll::PLL
    ωsys::Float64
    network_idx::UnitRange{Int64}
    filter_idx::UnitRange{Int64}
    pll_idx::UnitRange{Int64}
end

function run_pll_filter(network_file)
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

    # Convert filter capacitor voltage to RI for PLL initialization
    # This is the voltage that the PLL would measure
    Vd_flt = filter_states[VD_FLT]
    Vq_flt = filter_states[VQ_FLT]
    V_flt_ri = dq_ri(0.0) * [Vd_flt; Vq_flt]  # Convert DQ to RI

    # Initialize PLL states with filter capacitor voltage
    pll_states = initialize_pll(pll, V_flt_ri)

    # Combine all states
    states = vcat(network_states, filter_states, pll_states)

    # Define state indices
    network_idx = 1:NUM_STATES
    filter_idx = (NUM_STATES+1):(NUM_STATES+6)
    pll_idx = (NUM_STATES+6+1):(NUM_STATES+6+3)

    # Create parameter struct
    p = PLLFilterParams(
        network,
        filter,
        pll,
        ωsys,
        network_idx,
        filter_idx,
        pll_idx
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

    # Filter and PLL states are differential equations
    M_system[filter_idx] .= 1.0
    M_system[pll_idx] .= 1.0

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

    # System dynamics function
    function pll_filter_dynamics!(du, u, params, t)
        # Extract states for each component
        network_states = u[params.network_idx]
        filter_states = u[params.filter_idx]
        pll_states = u[params.pll_idx]

        # Make copies of states for Float64 compatibility
        network_states_f64 = convert.(Float64, network_states)
        filter_states_f64 = convert.(Float64, filter_states)
        pll_states_f64 = convert.(Float64, pll_states)

        # Arrays for derivatives
        du_network = similar(network_states_f64, length(params.network_idx))
        du_filter = similar(filter_states_f64, length(params.filter_idx))
        du_pll = similar(pll_states_f64, length(params.pll_idx))

        # Get terminal voltage and current from network at bus 2 (inverter bus)
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]
        v_grid = [v_2_d, v_2_q]
        
        # Get filter capacitor voltage for PLL
        v_flt_d = filter_states_f64[VD_FLT]
        v_flt_q = filter_states_f64[VQ_FLT]
        
        # Convert filter capacitor voltage to rectangular coordinates for PLL
        v_flt_ri = dq_ri(0.0) * [v_flt_d; v_flt_q]

        # Update PLL states
        update_pll_states!(
            pll_states_f64,
            du_pll,
            v_flt_ri,
            params.ωsys,
            params.pll
        )

        # Extract PLL angle
        θ_pll = pll_states_f64[THETA_IDX]
        
        # For this model, we'll just use a fixed inverter voltage, as we are waiting to implement the voltage controller (TODO: check this)
        v_inv = [Vd_inv, Vq_inv]  # Using initial values from filter initialization
        
        # Update filter states
        update_filter_states!(
            filter_states_f64,
            du_filter,
            v_inv,
            v_grid,
            params.ωsys,
            params.filter
        )

        # Get grid-side currents from filter
        i_2_d = filter_states_f64[ID_GRD]
        i_2_q = filter_states_f64[IQ_GRD]
        
        # Compute apparent power for the network
        P_terminal = v_2_d * i_2_d + v_2_q * i_2_q
        Q_terminal = v_2_q * i_2_d - v_2_d * i_2_q
        S_terminal = complex(P_terminal, Q_terminal)
        
        # Update network states
        _, _, _ = update_network_states!(
            network_states_f64,
            du_network,
            S_terminal,
            params.network
        )

        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.filter_idx] .= du_filter
        du[p.pll_idx] .= du_pll

        # Print some debugging info at integer time steps
        if abs(t - round(t)) < 0.00001
            println("t=$t: theta_pll=$(pll_states_f64[THETA_IDX]), vq_pll=$(pll_states_f64[VQ_PLL_IDX])")
        end
    end

    # Build function with mass matrix
    explicitDAE_M = ODEFunction(pll_filter_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 10.0)
    prob = ODEProblem(explicitDAE_M, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [5.0]

    # Define the condition for which to apply a perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end

    # Define the perturbation
    function affect!(integrator)
        # Choose one perturbation scenario
        
        # 1. Line Trip
        integrator.p.network.R_12 = 1e6
        integrator.p.network.X_12 = 1e6
        integrator.p.network.B_1 *= 0.5
        integrator.p.network.B_2 *= 0.5
        
        # 2. Frequency change
        # integrator.p.ωsys = 1.02  # 2% frequency increase
        
        # 3. Load change
        # integrator.p.network.Z_L *= 0.85  # 15% load increase
    end

    # Create a Callback function that represents the perturbation
    cb = DiscreteCallback(condition, affect!)

    # Run simulation
    sol = solve(prob, Rodas5(autodiff=false), dt=0.001, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    t = sol.t

    # Create plots
    
    # Network voltages
    p1 = plot(t, [sol[network_idx[V_2_D_IDX], i] for i in 1:length(t)],
        label="Vd", title="Bus 2 Voltage", linewidth=2)
    plot!(p1, t, [sol[network_idx[V_2_Q_IDX], i] for i in 1:length(t)],
        label="Vq", linewidth=2)
    
    # Filter currents
    p2 = plot(t, [sol[filter_idx[ID_INV], i] for i in 1:length(t)],
        label="Id_inv", title="Filter Currents", linewidth=2)
    plot!(p2, t, [sol[filter_idx[IQ_INV], i] for i in 1:length(t)],
        label="Iq_inv", linewidth=2)
    plot!(p2, t, [sol[filter_idx[ID_GRD], i] for i in 1:length(t)],
        label="Id_grd", linewidth=2)
    plot!(p2, t, [sol[filter_idx[IQ_GRD], i] for i in 1:length(t)],
        label="Iq_grd", linewidth=2)
    
    # Filter voltages
    p3 = plot(t, [sol[filter_idx[VD_FLT], i] for i in 1:length(t)],
        label="Vd_flt", title="Filter Capacitor Voltage", linewidth=2)
    plot!(p3, t, [sol[filter_idx[VQ_FLT], i] for i in 1:length(t)],
        label="Vq_flt", linewidth=2)
    
    # PLL states
    p4 = plot(t, [sol[pll_idx[THETA_IDX], i] for i in 1:length(t)],
        label="θ_pll", title="PLL Angle", linewidth=2)
    
    p5 = plot(t, [sol[pll_idx[VQ_PLL_IDX], i] for i in 1:length(t)],
        label="vq_pll", title="PLL q-axis Voltage", linewidth=2)
    
    # Line currents
    p6 = plot(t, [sol[network_idx[I_12_D_IDX], i] for i in 1:length(t)], 
        title="Line Currents", label="I_12_D", linewidth=2)
    plot!(p6, t, [sol[network_idx[I_12_Q_IDX], i] for i in 1:length(t)], 
        label="I_12_Q", linewidth=2)
    
    # Combine plots
    p_combined = plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(1200, 1800), left_margin=10mm)
    savefig(p_combined, "pll_filter_results.png")

    return sol
end

sol = run_pll_filter("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")