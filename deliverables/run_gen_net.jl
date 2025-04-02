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
const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6
const PSI_D = 7
const PSI_Q = 8

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

    network_states, _, _ = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)
    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)

    println("\nInitial States:")
    println("Machine states:")
    println("  Delta (rotor angle): $(machine_states[DELTA])")
    println("  Omega (rotor speed): $(machine_states[OMEGA])")
    println("  EQ_P: $(machine_states[EQ_P])")
    println("  ED_P: $(machine_states[ED_P])")
    println("  PSI_D_PP: $(machine_states[PSI_D_PP])")
    println("  PSI_Q_PP: $(machine_states[PSI_Q_PP])")
    println("  PSI_Q: $(machine_states[PSI_Q])")
    println("  PSI_D: $(machine_states[PSI_D])")

    println("Initial field voltage (Vf): $Vf_init")
    println("Initial mechanical torque (τm): $τ_m_init")

    states = vcat(network_states, machine_states)

    # Define state 
    network_idx = 1:22
    machine_idx = 23:30

    p = MachineNetworkParams(
        network,
        machine,
        Vf_init,
        τ_m_init,
        network_idx,
        machine_idx
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

    M_system[machine_idx] .= 1.0

    mass_matrix = Diagonal(M_system)
    println("Mass matrix = $(M_system)")

    # Define auxiliary variabls
    V_terminal_aux = Complex{Float64}[]
    S_terminal_aux = Complex{Float64}[]
    push!(V_terminal_aux, V_terminal)
    push!(S_terminal_aux, complex(P, Q))

    println("Initial Line Currents:")
    println("I_12_D: $(network_states[I_12_D_IDX])")
    println("I_12_Q: $(network_states[I_12_Q_IDX])")
    println("I_23_D: $(network_states[I_23_D_IDX])")
    println("I_23_Q: $(network_states[I_23_Q_IDX])")
    println("I_13_D: $(network_states[I_13_D_IDX])")
    println("I_13_Q: $(network_states[I_13_Q_IDX])")

    function machine_network_dynamics!(du, u, params, t)
        # Extract states for each component
        network_states = u[params.network_idx]
        machine_states = u[params.machine_idx]

        # Make copies of states for Float64 compatibility
        network_states_f64 = convert.(Float64, network_states)
        machine_states_f64 = convert.(Float64, machine_states)

        # Arrays for derivatives
        du_network = similar(network_states_f64, length(params.network_idx))
        du_machine = similar(machine_states_f64, length(params.machine_idx))

        # Grab terminal voltage from network
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]
        V_RI = dq_ri(0.0) * [v_2_d; v_2_q]
        V_terminal = complex(V_RI[1], V_RI[2])


        # Update the states of each component
        I_terminal_machine_pos, S_terminal_machine, ω_machine, V_mag, I_mag = update_machine_states!(
            machine_states_f64,
            du_machine,
            V_terminal,
            params.Vf,
            params.τm,
            params.machine
        )

        _, _, _ = update_network_states!(
            network_states_f64,
            du_network,
            S_terminal_machine,
            params.network
        )

        # Update auxiliary variables
        push!(V_terminal_aux, V_terminal)
        push!(S_terminal_aux, S_terminal_machine)

        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.machine_idx] .= du_machine

        if abs(t - round(t)) < 0.00001
            println("t=$t: δ=$(machine_states[DELTA]), ω=$ω_machine, τm=$(params.τm), Vf=$(params.Vf), V_mag=$V_mag, I_mag=$I_mag")
        end
    end

    # Build function 
    explicitDAE_M = ODEFunction(machine_network_dynamics!, mass_matrix=mass_matrix)

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
    sol = solve(prob, Rodas5(autodiff=false), dt=0.001, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

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

    p4 = plot(t, [sol[network_idx[I_12_D_IDX], i] for i in 1:length(t)], title="Line Currents", label="I_12_D", linewith=2)
    plot!(p4, t, [sol[network_idx[I_12_Q_IDX], i] for i in 1:length(t)], label="I_12_Q", linewith=2)
    plot!(p4, t, [sol[network_idx[I_23_D_IDX], i] for i in 1:length(t)], label="I_23_D", linewith=2)
    plot!(p4, t, [sol[network_idx[I_23_Q_IDX], i] for i in 1:length(t)], label="I_23_Q", linewith=2)
    plot!(p4, t, [sol[network_idx[I_13_D_IDX], i] for i in 1:length(t)], label="I_13_D", linewith=2)
    plot!(p4, t, [sol[network_idx[I_13_Q_IDX], i] for i in 1:length(t)], label="I_13_Q", linewith=2)

    p_combined = plot(p1, p2, p3, p4, layout=(4, 1), size=(1000, 1800))
    savefig(p_combined, "machine_network_results.png")

    return sol
end

sol = run_machine_network("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")