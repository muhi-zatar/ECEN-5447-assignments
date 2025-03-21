cd(@__DIR__)
using Pkg
Pkg.activate("../.")

using DifferentialEquations
using Plots
using LinearAlgebra
#using Sundials

include("NetworkModel.jl")
include("PowerFlowWrapper.jl")

using .NetworkModel
using .PowerFlowWrapper


mutable struct NetworkParams
    network::ThreeBusNetwork
    S_terminal::Complex{Float64}
    E_2_D::Float64
    E_2_Q::Float64
    network_idx::UnitRange{Int64}
    bus_idx::UnitRange{Int64}
end

function run_machine_network(network_file)
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
    network_idx = 1:16
    bus_idx = 17:18

    p = NetworkParams(
        network,
        complex(P, Q),
        0.0,
        0.0,
        network_idx,
        bus_idx
    )

    # Calculate EMF behind reactance at Bus 2
    p.E_2_D = network_states[9] - i_2_q_init * p.network.X_IB
    p.E_2_Q = network_states[10] - i_2_d_init * p.network.X_IB

    # M_system = zeros(Float64, length(states))

    # if isa(network.M, Vector)
    #     M_system[network_idx] .= network.M
    # else
    #     for i in 1:length(network_idx)
    #         M_system[network_idx[i]] = network.M[i, i]
    #     end
    # end

    # M_system[machine_idx] .= 1.0

    # mass_matrix = Diagonal(M_system)

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
        network_states = u[1:16]
        bus_2_states = u[17:18]

        network_states_f64 = convert.(Float64, network_states)
        bus_2_states_f64 = convert.(Float64, bus_2_states)

        # Arrays for derivatives
        du_network = similar(network_states, 16)
        du_bus_2 = similar(bus_2_states, 2)

        # Update the states
        V_terminal, S_terminal, I_terminal, i_2_d, i_2_q = update_network_states!(
            network_states_f64,
            du_network,
            S_terminal_aux[end],
            params.network
        )

        # Additional states to make Bus 2 an infinite bus
        du_bus_2[p.network_idx[1]] = ((network.E_IB_D - network_states[p.network_idx[9]]) + network.X_IB * i_2_q) / network.X_IB
        du_bus_2[p.network_idx[2]] = ((network.E_IB_Q - network_states[p.network_idx[10]]) - network.X_IB * i_2_d) / network.X_IB


        # Sanity check
        sanity_check(S_terminal, params.S_terminal, "Returned power", false)

        # Update auxiliary variables
        push!(t_aux, t)
        push!(V_terminal_aux, V_terminal)
        push!(S_terminal_aux, S_terminal)
        push!(I_terminal_aux, I_terminal)

        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.bus_idx] .= du_bus_2

        if abs(t - round(t)) < 0.00001
            println("t=$t: V_mag=$(abs(V_terminal)), I_mag=$(abs(Complex(i_2_d, i_2_q)))")
        end
    end

    # # Build function 
    # explicitDAE_M = ODEFunction(machine_network_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 10.0)
    prob = ODEProblem(network_dynamics!, states, tspan, p)

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
    sol = solve(prob, Tsit5(), dt=0.00005, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    p1 = plot(sol, idxs=(0, [7, 8, 9, 10, 11, 12]), title="Bus Voltages", labels=["V_1_D" "V_1_Q" "V_2_D" "V_2_Q" "V_3_D" "V_3_Q"])
    plot!(p1, t_aux, abs.(V_terminal_aux), label="Bus 2 Magnitude")
    plot!(twinx(), t_aux, rad2deg.(angle.(V_terminal_aux)), label="Bus 2 Angle", ylabel="Angle (degrees)", legend=:bottomleft, lw=2)
    p2 = plot(sol, idxs=(0, [13, 14, 15, 16]), title="Bus Current Injections", labels=["I_1_D" "I_1_Q" "I_3_D" "I_3_Q"])
    #plot!(p2, t_aux, i_2_d_aux, label="I_2_D")
    #plot!(p2, t_aux, i_2_q_aux, label="I_2_Q")
    plot!(p2, t_aux, abs.(I_terminal_aux), label="I_mag")
    plot!(twinx(), t_aux, rad2deg.(angle.(I_terminal_aux)), label="Bus 2 Angle", ylabel="Angle (degrees)", legend=:bottomleft, lw=2)
    p3 = plot(t_aux, real(S_terminal_aux), label="P", title="Bus Power Injection")
    plot!(p3, t_aux, imag(S_terminal_aux), label="Q")
    p_combined = plot(p1, p2, p3, layout=(3, 1), size=(1200, 2000))
    savefig(p_combined, "network_results.png")

    return sol
end

sol = run_machine_network("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")