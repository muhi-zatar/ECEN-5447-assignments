cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LinearAlgebra

include("NetworkModel.jl")
include("inverter_model/PLLModel.jl")
include("PowerFlowWrapper.jl")

using .NetworkModel
using .PLLModel
using .PowerFlowWrapper


# Definining some constants
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

# PLL state indices
const VQ_PLL_IDX = 1   # PLL q-axis voltage state (vq,pll)
const EPSILON_IDX = 2  # PLL error integral state (εpll)
const THETA_IDX = 3    # PLL angle (θpll)

# Matrix transformation functions
# These are used to convert between rectangular and polar coordinates
function ri_dq(δ::T) where {T<:Number}
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return Float64[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

function dq_ri(δ::T) where {T<:Number}
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return Float64[
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

mutable struct PllNetworkParams
    network::ThreeBusNetwork
    pll::PLL
    ωsys::Float64
    network_idx::UnitRange{Int64}
    pll_idx::UnitRange{Int64}
end

function run_pll_network(network_file)
    V_sol, θ_sol, P_sol, Q_sol = get_powerflow_results(network_file)

    V_mag = V_sol[BUS_MACHINE_MODEL]
    V_angle = θ_sol[BUS_MACHINE_MODEL]
    P = P_sol[BUS_MACHINE_MODEL]
    Q = Q_sol[BUS_MACHINE_MODEL]

    V_terminal = V_mag * exp(im * V_angle)
    println("Initial power flow results:")
    println("V_terminal = $V_mag p.u ∠ $V_angle rad")
    println("P = $P pu, Q = $Q pu")

    # Initialize the network
    network = ThreeBusNetwork()

    # Initialize the PLL
    pll = PLL()

    # Initialize network states
    network_states, Id_2, Iq_2 = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)

    # Extract the terminal voltage at the inverter bus (in network DQ) from the network states
    Vd_2 = network_states[V_2_D_IDX]
    Vq_2 = network_states[V_2_Q_IDX]

    # Convert DQ voltage to RI (rectangular) for PLL initialization
    V_ri = dq_ri(0.0) * [Vd_2; Vq_2]
    println("V_ri = $V_ri")
    # Initialize PLL states with grid voltage
    pll_states = initialize_pll(pll, V_ri)

    # Define system frequency (per unit)
    ωsys = 1.0

    # Combine all states
    states = vcat(network_states, pll_states)

    # Define state indices
    network_idx = 1:NUM_STATES
    pll_idx = (NUM_STATES+1):(NUM_STATES+3)

    # Create parameter struct
    p = PllNetworkParams(
        network,
        pll,
        ωsys,
        network_idx,
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

    # PLL states are differential equations
    M_system[pll_idx] .= 1.0

    mass_matrix = Diagonal(M_system)
    println("Mass matrix = $(M_system)")

    println("Initial Line Currents:")
    println("I_12_D: $(network_states[I_12_D_IDX])")
    println("I_12_Q: $(network_states[I_12_Q_IDX])")
    println("I_23_D: $(network_states[I_23_D_IDX])")
    println("I_23_Q: $(network_states[I_23_Q_IDX])")
    println("I_13_D: $(network_states[I_13_D_IDX])")
    println("I_13_Q: $(network_states[I_13_Q_IDX])")

    println("Initial PLL states:")
    println("vq_pll: $(pll_states[VQ_PLL_IDX])")
    println("epsilon_pll: $(pll_states[EPSILON_IDX])")
    println("theta_pll: $(pll_states[THETA_IDX])")

    function pll_network_dynamics!(du, u, params, t)
        # Extract states for each component
        network_states = u[params.network_idx]
        pll_states = u[params.pll_idx]

        # Make copies of states for Float64 compatibility
        network_states_f64 = convert.(Float64, network_states)
        pll_states_f64 = convert.(Float64, pll_states)

        # Arrays for derivatives
        du_network = similar(network_states_f64, length(params.network_idx))
        du_pll = similar(pll_states_f64, length(params.pll_idx))

        # Get terminal voltage from network at bus 2 (inverter bus)
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]

        # Convert network DQ voltage to rectangular coordinates for PLL
        v_2_ri = dq_ri(0.0) * [v_2_d; v_2_q]
        # println("v_2_ri = $v_2_ri")
        # Extract current PLL angle
        # θ_pll = pll_states_f64[THETA_IDX]
        
        # Update PLL states
        update_pll_states!(
            pll_states_f64,
            du_pll,
            v_2_ri,
            params.ωsys,
            params.pll
        )
        θ_pll = pll_states_f64[THETA_IDX]
        # convert back to network refernece framework to compute the power injection
        # Here we would typically have a current controller that uses the PLL angle so we need to use the filter
        # To make sure the problem is not in the reference frame, I will use at steady state the power flow solution.
        S_injection = complex(P, Q)  # Placeholder for power injection
        # println("S_injection = $S_injection")
        # exit(0)
        # Update network states with the power injection
        _, _, _ = update_network_states!(
            network_states_f64,
            du_network,
            S_injection,
            params.network
        )

        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.pll_idx] .= du_pll

        # Print some debugging info at integer time steps
        if abs(t - round(t)) < 0.00001
            println("t=$t: theta_pll=$(pll_states_f64[THETA_IDX]), vq_pll=$(pll_states_f64[VQ_PLL_IDX]), epsilon_pll=$(pll_states_f64[EPSILON_IDX])")
        end
    end

    # Build function with mass matrix
    explicitDAE_M = ODEFunction(pll_network_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 5.0)
    prob = ODEProblem(explicitDAE_M, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [25.0]

    # Define the condition for which to apply a perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end

    # Define the perturbation
    function affect!(integrator)
        #### Uncomment the desired perturbation ####
        # Load Jump
        # integrator.p.network.Z_L *= 0.85

        # Load Decrease
        # integrator.p.network.Z_L *= 1.15

        # Line Trip
        integrator.p.network.R_12 = 1e6
        integrator.p.network.X_12 = 1e6
        integrator.p.network.B_1 *= 0.5
        integrator.p.network.B_2 *= 0.5
        
        # Frequency change
        # integrator.p.ωsys = 1.02  # 2% frequency increase
    end

    # Create a Callback function that represents the perturbation
    cb = DiscreteCallback(condition, affect!)

    # Run simulation
    sol = solve(prob, Rodas5(autodiff=false), dt=0.001, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    t = sol.t

    # Plot PLL states
    p1 = plot(t, [sol[pll_idx[THETA_IDX], i] for i in 1:length(t)],
        label="θ_pll", title="PLL Angle", linewidth=2)
    
    p2 = plot(t, [sol[pll_idx[VQ_PLL_IDX], i] for i in 1:length(t)],
        label="vq_pll", title="PLL q-axis Voltage", linewidth=2)
    
    p3 = plot(t, [sol[pll_idx[EPSILON_IDX], i] for i in 1:length(t)],
        label="ε_pll", title="PLL Error Integral", linewidth=2)

    # Plot network voltages
    p4 = plot(t, [sol[network_idx[V_2_D_IDX], i] for i in 1:length(t)],
        label="Vd", title="Bus 2 Voltage", linewidth=2)
    plot!(p4, t, [sol[network_idx[V_2_Q_IDX], i] for i in 1:length(t)],
        label="Vq", linewidth=2)

    # Plot line currents
    p5 = plot(t, [sol[network_idx[I_12_D_IDX], i] for i in 1:length(t)], 
        title="Line Currents", label="I_12_D", linewidth=2)
    plot!(p5, t, [sol[network_idx[I_12_Q_IDX], i] for i in 1:length(t)], 
        label="I_12_Q", linewidth=2)
    plot!(p5, t, [sol[network_idx[I_23_D_IDX], i] for i in 1:length(t)], 
        label="I_23_D", linewidth=2)
    plot!(p5, t, [sol[network_idx[I_23_Q_IDX], i] for i in 1:length(t)], 
        label="I_23_Q", linewidth=2)

    p_combined = plot(p1, p2, p3, p4, p5, layout=(5, 1), size=(1000, 2000), left_margin=10mm)
    savefig(p_combined, "pll_network_results.png")

    return sol
end

sol = run_pll_network("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")