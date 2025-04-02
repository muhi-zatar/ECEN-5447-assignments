cd(@__DIR__)
using Pkg
Pkg.activate("../.")

using DifferentialEquations
using Plots
using LinearAlgebra

include("SauerPaiMachineModel.jl")

using .SauerPaiMachineModel

const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6
const PSI_D = 7
const PSI_Q = 8

mutable struct MachineParams
    machine::SauerPaiMachine
    V_terminal::Complex{Float64}  # Fixed terminal voltage
    Vf::Float64                  # Fixed field voltage
    τm::Float64                  # Fixed mechanical torque
end

function run_machine_only()
    V_mag = 1.0142
    V_angle = -0.007843738409384566
    P = 0.8
    Q = 0.10751144041475173

    V_terminal = V_mag * exp(im * V_angle)
    println("Initial conditions:")
    println("V_terminal = $V_mag p.u ∠ $V_angle rad")
    println("P = $P pu, Q = $Q pu")

    machine = SauerPaiMachine(
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)

    println("\nInitial Machine States:")
    println("Vf (field voltage): $Vf_init")
    println("Delta (rotor angle): $(machine_states[DELTA])")
    println("Omega (rotor speed): $(machine_states[OMEGA])")
    println("EQ_P: $(machine_states[EQ_P])")
    println("ED_P: $(machine_states[ED_P])")
    println("PSI_D_PP: $(machine_states[PSI_D_PP])")
    println("PSI_Q_PP: $(machine_states[PSI_Q_PP])")
    println("Initial field voltage (Vf): $Vf_init")
    println("Initial mechanical torque (τm): $τ_m_init")

    # Create params for ODE solver
    p = MachineParams(
        machine,
        V_terminal,
        Vf_init,
        τ_m_init
    )

    function machine_dynamics!(du, u, params, t)
        machine_states = u

        du_machine = zeros(Float64, 8)

        I_terminal_machine_pos, S_terminal_machine, ω_machine, V_mag, I_mag = update_machine_states!(
            machine_states,
            du_machine,
            params.V_terminal,
            params.Vf,
            params.τm,
            params.machine
        )

        du .= du_machine

        # if abs(t - round(t)) < 0.00001
        #     println("t=$t: δ=$(machine_states[DELTA]), ω=$ω_machine, τm=$(params.τm), Vf=$(params.Vf), V_mag=$V_mag, I_mag=$I_mag")
        # end
    end

    tspan = (0.0, 25.0)
    prob = ODEProblem(machine_dynamics!, machine_states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [2.5]               # Setting this far ahead for now – we can change this later

    # Define the condition for which to apply a perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end

    # Define the perturbation
    function affect!(integrator)
        #### Uncomment the desired perturbation ####
        # Load Jump
        integrator.p.τm *= 1.10

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

    p1 = plot(t, [sol[DELTA, i] for i in 1:length(t)],
        label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, [sol[OMEGA, i] for i in 1:length(t)],
        label="Rotor speed (ω)", linewidth=2)

    p2 = plot(t, [sol[EQ_P, i] for i in 1:length(t)],
        label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[ED_P, i] for i in 1:length(t)],
        label="ed'", linewidth=2)
    plot!(p2, t, [sol[PSI_D_PP, i] for i in 1:length(t)],
        label="ψd''", linewidth=2)
    plot!(p2, t, [sol[PSI_Q_PP, i] for i in 1:length(t)],
        label="ψq''", linewidth=2)

    p_combined = plot(p1, p2, layout=(2, 1), size=(800, 1200))
    savefig(p_combined, "machine_only_results.png")

    return sol
end

sol = run_machine_only()

println("Simulation complete!")