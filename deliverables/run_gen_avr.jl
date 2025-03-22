cd(@__DIR__)
using Pkg
Pkg.activate("../.")

using DifferentialEquations
using Plots
using LinearAlgebra

include("SauerPaiMachineModel.jl")
include("EXST1Model.jl")

using .SauerPaiMachineModel
using .EXST1Model

const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6

const VF_IDX = 1
const VT_IDX = 2
const VLL_IDX = 3
const VFB_IDX = 4

mutable struct MachineAVRParams
    machine::SauerPaiMachine
    avr::EXST1
    V_terminal::Complex{Float64}
    τm::Float64
    machine_idx::UnitRange{Int64}
    avr_idx::UnitRange{Int64}
end

function run_machine_avr()
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

    avr = EXST1(
        V_ref=V_mag  # Set reference voltage to initial voltage
    )

    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)
    avr_states = initialize_avr(avr, V_mag, Vf_init)

    println("\nInitial States:")
    println("Machine states:")
    println("  Delta (rotor angle): $(machine_states[DELTA])")
    println("  Omega (rotor speed): $(machine_states[OMEGA])")
    println("  EQ_P: $(machine_states[EQ_P])")
    println("  ED_P: $(machine_states[ED_P])")
    println("  PSI_D_PP: $(machine_states[PSI_D_PP])")
    println("  PSI_Q_PP: $(machine_states[PSI_Q_PP])")

    println("AVR states:")
    println("  Field Voltage: $(avr_states[VF_IDX])")
    println("  Terminal Voltage: $(avr_states[VT_IDX])")
    println("  Lead-Lag: $(avr_states[VLL_IDX])")
    println("  Feedback: $(avr_states[VFB_IDX])")

    println("Initial field voltage (Vf): $Vf_init")
    println("Initial mechanical torque (τm): $τ_m_init")

    states = vcat(machine_states, avr_states)

    machine_idx = 1:6
    avr_idx = 7:10

    p = MachineAVRParams(
        machine,
        avr,
        V_terminal,
        τ_m_init,
        machine_idx,
        avr_idx
    )

    # Define mass matrix for system
    M_system = zeros(Float64, length(states))

    M_system[machine_idx] .= machine.M
    M_system[avr_idx] .= 1.0

    println("\nDerivative Coefficients: $M_system\n")

    mass_matrix = Diagonal(M_system)


    function machine_avr_dynamics!(du, u, params, t)
        machine_states = @view u[params.machine_idx]
        avr_states = @view u[params.avr_idx]

        du_machine = zeros(Float64, length(params.machine_idx))
        du_avr = zeros(Float64, length(params.avr_idx))

        Vf = avr_states[VF_IDX]

        I_terminal_machine_pos, S_terminal_machine, ω_machine, V_mag, I_mag = update_machine_states!(
            machine_states,
            du_machine,
            params.V_terminal,
            Vf,
            params.τm,
            params.machine
        )

        update_avr_states!(
            avr_states,
            du_avr,
            V_mag,
            params.avr
        )

        du[params.machine_idx] .= du_machine
        du[params.avr_idx] .= du_avr

        if abs(t - round(t)) < 0.001
            println("t=$t: δ=$(machine_states[DELTA]), ω=$ω_machine, τm=$(params.τm), Vf=$(Vf), V_mag=$V_mag, I_mag=$I_mag")
        end
    end

    f = ODEFunction(machine_avr_dynamics!, mass_matrix=mass_matrix)

    tspan = (0.0, 20.0)
    prob = ODEProblem(f, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [25.0]               # Setting this far ahead for now – we can change this later

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
    sol = solve(prob, Rodas5P(autodiff=false), saveat=0.01, callback=cb, tstops=perturb_times)

    t = sol.t

    p1 = plot(t, [sol[machine_idx[DELTA], i] for i in 1:length(t)],
        label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, [sol[machine_idx[OMEGA], i] for i in 1:length(t)],
        label="Rotor speed (ω)", linewidth=2)

    p2 = plot(t, [sol[machine_idx[EQ_P], i] for i in 1:length(t)],
        label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[machine_idx[ED_P], i] for i in 1:length(t)],
        label="ed'", linewidth=2)
    plot!(p2, t, [sol[machine_idx[PSI_D_PP], i] for i in 1:length(t)],
        label="ψd''", linewidth=2)
    plot!(p2, t, [sol[machine_idx[PSI_Q_PP], i] for i in 1:length(t)],
        label="ψq''", linewidth=2)

    p3 = plot(t, [sol[avr_idx[VF_IDX], i] for i in 1:length(t)],
        label="Field Voltage (Vf)", title="AVR States", linewidth=2)
    plot!(p3, t, [sol[avr_idx[VT_IDX], i] for i in 1:length(t)],
        label="Terminal Voltage (Vt)", linewidth=2)
    plot!(p3, t, [sol[avr_idx[VLL_IDX], i] for i in 1:length(t)],
        label="Lead-Lag State", linewidth=2)
    plot!(p3, t, [sol[avr_idx[VFB_IDX], i] for i in 1:length(t)],
        label="Feedback State", linewidth=2)

    τe_values = []
    for i in 1:length(t)
        # Simplified approximation of electrical torque
        delta = sol[machine_idx[DELTA], i]
        eq_p = sol[machine_idx[EQ_P], i]
        ed_p = sol[machine_idx[ED_P], i]

        # This is a simplified approximation
        τe = eq_p * sin(delta) - ed_p * cos(delta)
        push!(τe_values, τe)
    end

    p4 = plot(t, τe_values,
        label="Electrical Torque (approx)", title="Machine Torques", linewidth=2)
    plot!(p4, t, fill(p.τm, length(t)),
        label="Mechanical Torque (constant)", linewidth=2)

    p_combined = plot(p1, p2, p3, p4, layout=(4, 1), size=(800, 1600))
    savefig(p_combined, "machine_avr_results_new.png")

    return sol
end

sol = run_machine_avr()

println("Simulation complete!")