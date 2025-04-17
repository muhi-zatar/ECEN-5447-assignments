cd(@__DIR__)
using Pkg
Pkg.activate(".")

using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LinearAlgebra

include("NetworkModel.jl")
include("inverter_model/FilterModel.jl")
include("PowerFlowWrapper.jl")

using .NetworkModel
using .FilterModel
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

# Filter state indices
const NUM_FLT_STATES = 6
const ID_INV = 1       # Inverter-side inductor current (real component)
const IQ_INV = 2       # Inverter-side inductor current (imag component)
const VD_FLT = 3       # Filter capacitor voltage (real component)
const VQ_FLT = 4       # Filter capacitor voltage (imag component)
const ID_GRD = 5       # Grid-side inductor current (real component)
const IQ_GRD = 6       # Grid-side inductor current (imag component)

# Matrix transformation functions
# These are used to convert between rectangular and polar coordinates
# Same as power system package to avoid confusion in this part.

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

mutable struct FilterNetworkParams
    network::ThreeBusNetwork
    flt::Filter
    Vd_inv::Float64
    Vq_inv::Float64
    ωg::Float64
    network_idx::UnitRange{Int64}
    filter_idx::UnitRange{Int64}
end

function run_filter_network(network_file)
    V_sol, θ_sol, P_sol, Q_sol = get_powerflow_results(network_file)

    V_mag = V_sol[BUS_MACHINE_MODEL]
    V_angle = θ_sol[BUS_MACHINE_MODEL]
    P = P_sol[BUS_MACHINE_MODEL]
    Q = Q_sol[BUS_MACHINE_MODEL]

    V_terminal = V_mag * exp(im * V_angle)
    println("Initial power flow results:")
    println("V_terminal = $V_mag p.u ∠ $V_angle rad")
    println("P = $P pu, Q = $Q pu")

    network = ThreeBusNetwork()

    flt = Filter()

    network_states, Id_grd, Iq_grd = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)

    # Extract the terminal voltage at the inverter bus (in network DQ) from the network states
    Vd_grd = network_states[V_2_D_IDX]
    Vq_grd = network_states[V_2_Q_IDX]

    # Pack the terminal voltage and current (in network DQ) into a vector to pass to the filter
    v_term = [Vd_grd, Vq_grd]
    i_term = [Id_grd, Iq_grd]

    # Define a system frequency (I think this is meant to be per-unit based on Darco paper)
    ωg = 1.0

    filter_states, Vd_inv, Vq_inv = initialize_filter(flt, v_term, i_term, ωg)

    #println("\nInitial States:")
    #println("Filter states:")

    states = vcat(network_states, filter_states)

    # Define state 
    network_idx = 1:22
    filter_idx = 23:28

    p = FilterNetworkParams(
        network,
        flt,
        Vd_inv,
        Vq_inv,
        ωg,
        network_idx,
        filter_idx
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

    M_system[filter_idx] .= 1.0

    mass_matrix = Diagonal(M_system)
    println("Mass matrix = $(M_system)")

    println("Initial Line Currents:")
    println("I_12_D: $(network_states[I_12_D_IDX])")
    println("I_12_Q: $(network_states[I_12_Q_IDX])")
    println("I_23_D: $(network_states[I_23_D_IDX])")
    println("I_23_Q: $(network_states[I_23_Q_IDX])")
    println("I_13_D: $(network_states[I_13_D_IDX])")
    println("I_13_Q: $(network_states[I_13_Q_IDX])")

    function filter_network_dynamics!(du, u, params, t)
        # Extract states for each component
        network_states = u[params.network_idx]
        filter_states = u[params.filter_idx]

        # Make copies of states for Float64 compatibility
        network_states_f64 = convert.(Float64, network_states)
        filter_states_f64 = convert.(Float64, filter_states)

        # Arrays for derivatives
        du_network = similar(network_states_f64, length(params.network_idx))
        du_filter = similar(filter_states_f64, length(params.filter_idx))

        # Grab terminal voltage from network
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]
        v_grid = [v_2_d, v_2_q]

        # Pack a vector of the inverter voltage
        v_inv = [params.Vd_inv, params.Vq_inv]

        # Grab the inverter current from the filter
        i_2_d = filter_states_f64[ID_GRD]
        i_2_q = filter_states_f64[IQ_GRD]

        # Compute apparent power for the network to use
        P_terminal = v_2_d * i_2_d + v_2_q * i_2_q
        Q_terminal = v_2_q * i_2_d - v_2_d * i_2_q
        S_terminal = complex(P_terminal, Q_terminal)

        # Update the states of each component
        update_filter_states!(filter_states_f64, du_filter, v_inv, v_grid, params.ωg, params.flt)

        _, _, _ = update_network_states!(
            network_states_f64,
            du_network,
            S_terminal,
            params.network
        )

        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.filter_idx] .= du_filter

        if abs(t - round(t)) < 0.00001
            println("t=$t: Id_inv=$(filter_states_f64[ID_INV]), Id_grd=$(filter_states_f64[ID_GRD]), S_terminal=$S_terminal")
        end
    end

    # Build function 
    explicitDAE_M = ODEFunction(filter_network_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 20.0)
    prob = ODEProblem(explicitDAE_M, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [2.5]

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
    end

    # Create a Callback function that represents the perturbation
    cb = DiscreteCallback(condition, affect!)

    # Run simulation
    sol = solve(prob, Rodas5(autodiff=false), dt=0.001, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    t = sol.t

    p1 = plot(t, [sol[filter_idx[ID_INV], i] for i in 1:length(t)],
        label="Id (inv)", title="Filter Currents", linewidth=2)
    plot!(p1, t, [sol[filter_idx[IQ_INV], i] for i in 1:length(t)],
        label="Iq (inv)", linewidth=2)
    plot!(p1, t, [sol[filter_idx[ID_GRD], i] for i in 1:length(t)],
        label="Iq (grid)", linewidth=2)
    plot!(p1, t, [sol[filter_idx[IQ_GRD], i] for i in 1:length(t)],
        label="Iq (grid)", linewidth=2)

    p2 = plot(t, [sol[filter_idx[VD_FLT], i] for i in 1:length(t)],
        label="Vd", title="Filter Capacitor Voltages", linewidth=2)
    plot!(p2, t, [sol[filter_idx[VQ_FLT], i] for i in 1:length(t)],
        label="Vq", linewidth=2)

    p3 = plot(t, [sol[network_idx[I_12_D_IDX], i] for i in 1:length(t)], title="Line Currents", label="I_12_D", linewith=2)
    plot!(p3, t, [sol[network_idx[I_12_Q_IDX], i] for i in 1:length(t)], label="I_12_Q", linewith=2)
    plot!(p3, t, [sol[network_idx[I_23_D_IDX], i] for i in 1:length(t)], label="I_23_D", linewith=2)
    plot!(p3, t, [sol[network_idx[I_23_Q_IDX], i] for i in 1:length(t)], label="I_23_Q", linewith=2)
    plot!(p3, t, [sol[network_idx[I_13_D_IDX], i] for i in 1:length(t)], label="I_13_D", linewith=2)
    plot!(p3, t, [sol[network_idx[I_13_Q_IDX], i] for i in 1:length(t)], label="I_13_Q", linewith=2)

    p_combined = plot(p1, p2, p3, layout=(3, 1), size=(1000, 1800), left_margin=10mm)
    savefig(p_combined, "filter_network_results.png")

    return sol
end

sol = run_filter_network("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")