cd(@__DIR__)
using Pkg
Pkg.activate("../.")

using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LinearAlgebra
#using Sundials

include("NetworkModel.jl")
include("PowerFlowWrapper.jl")

using .NetworkModel
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


mutable struct NetworkParams
    network::ThreeBusNetwork
    S_terminal::Complex{Float64}
    E_2_D::Float64
    E_2_Q::Float64
    network_idx::UnitRange{Int64}
    bus_idx::UnitRange{Int64}
end

function run_network(network_file)
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

    network_states, i_2_d_init, i_2_q_init = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)
    I_terminal_init = (i_2_d_init + im * i_2_q_init) * exp(-im * π / 2)

    states = vcat(network_states, i_2_d_init, i_2_q_init)
    network_idx = 1:22
    bus_idx = 23:24

    p = NetworkParams(
        network,
        complex(P, Q),
        0.0,
        0.0,
        network_idx,
        bus_idx
    )

    # Calculate EMF behind reactance at Bus 2
    p.E_2_D = network_states[V_2_D_IDX] - i_2_q_init * p.network.X_IB
    p.E_2_Q = network_states[V_2_Q_IDX] + i_2_d_init * p.network.X_IB

    # Build mass matrix
    M_system = zeros(Float64, length(states))

    if isa(network.M, Vector)
        M_system[network_idx] .= network.M
    else
        for i in 1:length(network_idx)
            M_system[network_idx[i]] = network.M[i, i]
        end
    end

    M_system[bus_idx] .= 0.0            # Infinite bus current injection is algebraic

    mass_matrix = Diagonal(M_system)

    println("Mass Matrix: $(M_system)")

    # Define auxiliary variabls
    t_aux = Float64[]
    V_terminal_aux = Complex{Float64}[]
    S_terminal_aux = Complex{Float64}[]
    I_terminal_aux = Complex{Float64}[]
    i_1_d_aux = Float64[]
    i_1_q_aux = Float64[]
    i_2_d_aux = Float64[]
    i_2_q_aux = Float64[]
    i_3_d_aux = Float64[]
    i_3_q_aux = Float64[]
    push!(t_aux, 0.0)
    push!(V_terminal_aux, V_terminal)
    push!(S_terminal_aux, complex(P, Q))
    push!(I_terminal_aux, I_terminal_init)
    push!(i_1_d_aux, network_states[I_12_D_IDX])
    push!(i_1_q_aux, network_states[I_12_Q_IDX])
    push!(i_2_d_aux, network_states[I_23_D_IDX])
    push!(i_2_q_aux, network_states[I_23_D_IDX])
    push!(i_3_d_aux, network_states[I_13_D_IDX])
    push!(i_3_q_aux, network_states[I_13_Q_IDX])

    function network_dynamics!(du, u, params, t)
        network_states = u[params.network_idx]
        bus_2_states = u[params.bus_idx]

        network_states_f64 = convert.(Float64, network_states)
        bus_2_states_f64 = convert.(Float64, bus_2_states)

        # Arrays for derivatives
        du_network = similar(network_states, length(network_states))
        du_bus_2 = similar(bus_2_states, length(bus_2_states))

        # Calculate new complex power
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]
        i_2_d = bus_2_states_f64[1]
        i_2_q = bus_2_states_f64[2]
        P_new = v_2_d * i_2_d + v_2_q * i_2_q
        Q_new = v_2_q * i_2_d - v_2_d * i_2_q
        S_new = complex(P_new, Q_new)

        # Update the states
        # Note that the returned values here are only used for plotting via the auxiliary variables
        V_terminal, _, I_terminal, _, _ = update_network_states!(
            network_states_f64,
            du_network,
            S_new,
            params.network
        )

        # Additional states to make Bus 2 an infinite bus
        du_bus_2[p.network_idx[1]] = ((params.E_2_D - network_states[p.network_idx[9]]) + network.X_IB * i_2_q)
        du_bus_2[p.network_idx[2]] = ((params.E_2_Q - network_states[p.network_idx[10]]) - network.X_IB * i_2_d)

        # Update auxiliary variables
        push!(t_aux, t)
        push!(V_terminal_aux, V_terminal)
        push!(S_terminal_aux, S_new)
        push!(I_terminal_aux, I_terminal)
        push!(i_1_d_aux, network_states_f64[I_12_D_IDX])
        push!(i_1_q_aux, network_states_f64[I_12_Q_IDX])
        push!(i_2_d_aux, network_states_f64[I_23_D_IDX])
        push!(i_2_q_aux, network_states_f64[I_23_Q_IDX])
        push!(i_3_d_aux, network_states_f64[I_13_D_IDX])
        push!(i_3_q_aux, network_states_f64[I_13_Q_IDX])

        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.bus_idx] .= du_bus_2

        if abs(t - round(t)) < 0.00001
            println("t=$t: V_mag=$(abs(V_terminal)), I_mag=$(abs(Complex(i_2_d, i_2_q)))")
        end
    end

    # Build function 
    explicitDAE_M = ODEFunction(network_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 30.0)
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
        #integrator.p.network.Z_L *= 0.85

        # Load Decrease
        #integrator.p.network.Z_L *= 1.15

        # Line Trip
        integrator.p.network.R_12 = 1e6
        integrator.p.network.X_12 = 1e6
        integrator.p.network.B_1 *= 0.5
        integrator.p.network.B_2 *= 0.5
    end

    # Create a Callback function that represents the perturbation
    cb = DiscreteCallback(condition, affect!)

    # Run simulation
    sol = solve(prob, Rosenbrock23(autodiff=false), dt=0.001, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    p1 = plot(sol, idxs=(0, [7, 8, 9, 10, 11, 12]), title="Bus Voltages", labels=["V_1_D" "V_1_Q" "V_2_D" "V_2_Q" "V_3_D" "V_3_Q"])
    plot!(p1, t_aux, abs.(V_terminal_aux), label="Bus 2 Magnitude")
    p2 = plot(t_aux, abs.(I_terminal_aux), label="I_mag", title="Bus Currents")
    plot!(p2, t_aux, i_1_d_aux, label="I_12_D")
    plot!(p2, t_aux, i_1_q_aux, label="I_12_Q")
    plot!(p2, t_aux, i_2_d_aux, label="I_23_D")
    plot!(p2, t_aux, i_2_q_aux, label="I_23_Q")
    plot!(p2, t_aux, i_3_d_aux, label="I_13_D")
    plot!(p2, t_aux, i_3_q_aux, label="I_13_Q")
    p3 = plot(t_aux, real(S_terminal_aux), label="P", title="Bus Power Injection")
    plot!(p3, t_aux, imag(S_terminal_aux), label="Q")
    p_combined = plot(p1, p2, p3, layout=(3, 1), size=(1200, 2000), left_margin=50mm)
    savefig(p_combined, "network_results_DAE.png")

    return sol
end

sol = run_network("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")