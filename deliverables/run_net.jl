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
const NUM_STATES = 16
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
# const I_1_D_IDX = 13
# const I_1_Q_IDX = 14
# const I_3_D_IDX = 15
# const I_3_Q_IDX = 16


mutable struct NetworkParams
    network::ThreeBusNetwork
    S_terminal::Complex{Float64}
    E_2_D::Float64
    E_2_Q::Float64
    network_idx::UnitRange{Int64}
end

function verify_network_initialiation(states, i_2_d, i_2_q, pf_solution)
    i_12_d = states[I_12_D_IDX]
    i_12_q = states[I_12_Q_IDX]
    i_13_d = states[I_13_D_IDX]
    i_13_q = states[I_13_Q_IDX]
    i_23_d = states[I_23_D_IDX]
    i_23_q = states[I_23_Q_IDX]
    v_1_d = states[V_1_D_IDX]
    v_1_q = states[V_1_Q_IDX]
    v_2_d = states[V_2_D_IDX]
    v_2_q = states[V_2_Q_IDX]
    v_3_d = states[V_3_D_IDX]
    v_3_q = states[V_3_Q_IDX]
    i_1_d = states[I_1_D_IDX]
    i_1_q = states[I_1_Q_IDX]
    i_3_d = states[I_3_D_IDX]
    i_3_q = states[I_3_Q_IDX]

    V_sol, θ_sol, P_sol, Q_sol = pf_solution

    P_1 = (v_1_d .* i_1_d) .+ (v_1_q .* i_1_q)
    Q_1 = (v_1_q .* i_1_d) .- (v_1_d .* i_1_q)
    P_2 = (v_2_d .* i_2_d) .+ (v_2_q .* i_2_q)
    Q_2 = (v_2_q .* i_2_d) .- (v_2_d .* i_2_q)
    P_3 = (v_3_d .* i_3_d) .+ (v_3_q .* i_3_q)
    Q_3 = (v_3_q .* i_3_d) .- (v_3_d .* i_3_q)

    P_test = [P_1; P_2; P_3]
    Q_test = [Q_1; Q_2; Q_3]

    @show P_test
    @show Q_test

    sanity_check(P_test, P_sol, "P")
    sanity_check(Q_test, Q_sol, "Q")

    return
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

    #### Sanity Check on Network Initialization ####
    #verify_network_initialiation(network_states, i_2_d_init, i_2_q_init, (V_sol, θ_sol, P_sol, Q_sol))

    states = vcat(network_states, i_2_d_init, i_2_q_init)
    network_idx = 1:12

    p = NetworkParams(
        network,
        complex(P, Q),
        0.0,
        0.0,
        network_idx
    )

    # Calculate EMF behind reactance at Bus 2
    p.E_2_D = network_states[V_2_D_IDX] - i_2_q_init * p.network.X_IB
    p.E_2_Q = network_states[V_2_Q_IDX] + i_2_d_init * p.network.X_IB

    # Define auxiliary variabls
    t_aux = Float64[]
    V_terminal_aux = Complex{Float64}[]
    S_terminal_aux = Complex{Float64}[]
    I_terminal_aux = Complex{Float64}[]
    push!(t_aux, 0.0)
    push!(V_terminal_aux, V_terminal)
    push!(S_terminal_aux, complex(P, Q))
    push!(I_terminal_aux, I_terminal_init)

    function network_dynamics!(du, u, params, t)
        network_states = u

        network_states_f64 = convert.(Float64, network_states)

        # Arrays for derivatives
        du_network = similar(network_states, 12)

        # Calculate new complex power
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]
        #println("V_dq @ 2: $(complex(v_2_d, v_2_q))")
        i_2_q = -1 * (params.E_2_D - v_2_d) / network.X_IB
        i_2_d = (params.E_2_Q - v_2_q) / network.X_IB
        #println("I_dq @ 2: $(complex(i_2_d, i_2_q))")
        P_new = v_2_d * i_2_d + v_2_q * i_2_q
        Q_new = v_2_q * i_2_d - v_2_d * i_2_q
        S_new = complex(P_new, Q_new)

        # Update the states
        V_terminal, S_terminal, I_terminal, i_2_d, i_2_q = update_network_states!(
            network_states_f64,
            du_network,
            S_new,
            params.network
        )

        # Sanity check
        #println("S_new = $S_new")
        #sanity_check(S_new, params.S_terminal, "Returned power", false)

        # Update auxiliary variables
        push!(t_aux, t)
        push!(V_terminal_aux, V_terminal)
        push!(S_terminal_aux, S_new)
        push!(I_terminal_aux, I_terminal)

        # Copy derivatives to output
        du[p.network_idx] .= du_network

        if abs(t - round(t)) < 0.00001
            println("t=$t: V_mag=$(abs(V_terminal)), I_mag=$(abs(Complex(i_2_d, i_2_q)))")
        end
    end

    # # Build function 
    # explicitDAE_M = ODEFunction(machine_network_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 20.0)
    prob = ODEProblem(network_dynamics!, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [5.0]

    # Define the condition for which to apply a perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end

    # Define the perturbation
    function affect!(integrator)
        #### Uncomment the desired perturbation ####
        # Load Jump
        integrator.p.network.Z_L *= 1.05

        # Load Decrease
        #integrator.p.network.Z_L *= 1.15

        # Line Trip
        #integrator.p.network.R_12 = 1e6
        #integrator.p.network.X_12 = 1e6
    end

    # Create a Callback function that represents the perturbation
    cb = DiscreteCallback(condition, affect!)

    # Run simulation
    sol = solve(prob, Tsit5(), dt=0.005, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    p1 = plot(sol, idxs=(0, [7, 8, 9, 10, 11, 12]), title="Bus Voltages", labels=["V_1_D" "V_1_Q" "V_2_D" "V_2_Q" "V_3_D" "V_3_Q"])
    plot!(p1, t_aux, abs.(V_terminal_aux), label="Bus 2 Magnitude")
    #plot!(twinx(), t_aux, rad2deg.(angle.(V_terminal_aux)), label="Bus 2 Angle", ylabel="Angle (degrees)", legend=:bottomleft, lw=2)
    #p2 = plot(sol, idxs=(0, [13, 14, 15, 16]), title="Bus Current Injections", labels=["I_1_D" "I_1_Q" "I_3_D" "I_3_Q"])
    #plot!(p2, t_aux, i_2_d_aux, label="I_2_D")
    #plot!(p2, t_aux, i_2_q_aux, label="I_2_Q")
    p2 = plot(t_aux, abs.(I_terminal_aux), label="I_mag")
    #plot!(twinx(), t_aux, rad2deg.(angle.(I_terminal_aux)), label="Bus 2 Angle", ylabel="Angle (degrees)", legend=:bottomleft, lw=2)
    p3 = plot(t_aux, real(S_terminal_aux), label="P", title="Bus Power Injection")
    plot!(p3, t_aux, imag(S_terminal_aux), label="Q")
    p_combined = plot(p1, p2, p3, layout=(3, 1), size=(1200, 2000), left_margin=50mm)
    savefig(p_combined, "network_results.png")

    return sol
end

sol = run_network("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")