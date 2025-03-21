cd(@__DIR__)
using Pkg
Pkg.activate("../.")

using DifferentialEquations
using Plots
using LinearAlgebra
#using Sundials

include("NetworkModel.jl")
include("SauerPaiMachineModel.jl")
include("PowerFlowWrapper.jl")

using .NetworkModel
using .SauerPaiMachineModel
using .PowerFlowWrapper


const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6

mutable struct MachineNetworkParams
    network::ThreeBusNetwork
    machine::SauerPaiMachine
    Vf::Float64
    τm::Float64
    network_idx::UnitRange{Int64}
    machine_idx::UnitRange{Int64}
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

    machine = SauerPaiMachine(
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

    network_states, i_2_d_init, i_2_q_init = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)
    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)

    println("\nInitial States:")
    println("Machine states:")
    println("  Delta (rotor angle): $(machine_states[DELTA])")
    println("  Omega (rotor speed): $(machine_states[OMEGA])")
    println("  EQ_P: $(machine_states[EQ_P])")
    println("  ED_P: $(machine_states[ED_P])")
    println("  PSI_D_PP: $(machine_states[PSI_D_PP])")
    println("  PSI_Q_PP: $(machine_states[PSI_Q_PP])")

    println("Initial field voltage (Vf): $Vf_init")
    println("Initial mechanical torque (τm): $τ_m_init")

    states = vcat(network_states, machine_states)

    # Define state 
    network_idx = 1:16
    machine_idx = 17:22

    p = MachineNetworkParams(
        network,
        machine,
        Vf_init,
        τ_m_init,
        network_idx,
        machine_idx
    )

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
    V_terminal_aux = Complex{Float64}[]
    S_terminal_aux = Complex{Float64}[]
    t_aux = Float64[]
    i_2_d_aux = Float64[]
    i_2_q_aux = Float64[]
    push!(V_terminal_aux, V_terminal)
    push!(S_terminal_aux, complex(P, Q))
    push!(t_aux, 0.0)
    push!(i_2_d_aux, i_2_d_init)
    push!(i_2_q_aux, i_2_q_init)

    function machine_network_dynamics!(du, u, params, t)
        # Extract states for each component
        network_states = u[params.network_idx]
        machine_states = u[params.machine_idx]

        # Make copies of states for Float64 compatibility
        network_states_f64 = convert.(Float64, network_states)
        machine_states_f64 = convert.(Float64, machine_states)

        # Arrays for derivatives
        du_network = similar(network_states, length(params.network_idx))
        du_machine = similar(machine_states, length(params.machine_idx))

        # Update the states of each component
        I_terminal_machine_pos, S_terminal_machine, ω_machine, V_mag, I_mag = update_machine_states!(
            machine_states_f64,
            du_machine,
            V_terminal_aux[end],
            params.Vf,
            params.τm,
            params.machine
        )

        V_terminal, S_terminal_network, I_terminal_network, i_2_d, i_2_q = update_network_states!(
            network_states_f64,
            du_network,
            S_terminal_aux[end],
            params.network
        )

        #sanity_check(S_terminal_network, S_terminal_machine, "Machine vs. Network Power", false)
        #sanity_check(I_terminal_network, I_terminal_machine_pos, "Machine vs. Network Current", false)

        # Update auxiliary variables
        push!(V_terminal_aux, V_terminal)
        push!(S_terminal_aux, S_terminal_machine)
        push!(t_aux, t)
        push!(i_2_d_aux, i_2_d)
        push!(i_2_q_aux, i_2_q)

        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.machine_idx] .= du_machine

        if abs(t - round(t)) < 0.001
            println("t=$t: δ=$(machine_states[DELTA]), ω=$ω_machine, τm=$(params.τm), Vf=$(params.Vf), V_mag=$V_mag, I_mag=$I_mag")
        end
    end

    # # Build function 
    # explicitDAE_M = ODEFunction(machine_network_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 20.0)
    prob = ODEProblem(machine_network_dynamics!, states, tspan, p)

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
    sol = solve(prob, Tsit5(), dt=0.00005, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    t = sol.t

    p1 = plot(t, [sol[machine_idx[1], i] for i in 1:length(t)],
        label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, [sol[machine_idx[2], i] for i in 1:length(t)],
        label="Rotor speed (ω)", linewidth=2)

    p2 = plot(t, [sol[machine_idx[3], i] for i in 1:length(t)],
        label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[machine_idx[4], i] for i in 1:length(t)],
        label="ed'", linewidth=2)
    plot!(p2, t, [sol[machine_idx[5], i] for i in 1:length(t)],
        label="ψd''", linewidth=2)
    plot!(p2, t, [sol[machine_idx[6], i] for i in 1:length(t)],
        label="ψq''", linewidth=2)

    p3 = plot(t, [sol[network_idx[9], i] for i in 1:length(t)],
        title="Terminal Voltage (Machine Bus)", label="Vd", linewidth=2)
    plot!(p3, t, [sol[network_idx[10], i] for i in 1:length(t)],
        label="Vq", linewidth=2)
    p4 = plot(sol, idxs=(0, [7, 8, 9, 10, 11, 12]), title="Network Bus Voltages", labels=["V_1_D" "V_1_Q" "V_2_D" "V_2_Q" "V_3_D" "V_3_Q"])
    p5 = plot(sol, idxs=(0, [13, 14, 15, 16]), title="Network Bus Current Injections", labels=["I_1_D" "I_1_Q" "I_3_D" "I_3_Q"])
    plot!(p5, t_aux, i_2_d_aux, label="I_2_D")
    plot!(p5, t_aux, i_2_q_aux, label="I_2_Q")
    p_combined = plot(p1, p2, p3, p4, p5, layout=(5, 1), size=(1200, 2000))
    savefig(p_combined, "machine_network_results.png")

    return sol
end

sol = run_machine_network("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")