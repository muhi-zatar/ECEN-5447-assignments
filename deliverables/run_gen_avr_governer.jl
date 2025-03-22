# Ensure script runs in the correct environment
cd(@__DIR__)
using Pkg
Pkg.activate("../.")

# Importing the required libraries
using DifferentialEquations
using Plots
using Plots.PlotMeasures
using LinearAlgebra

# Import machine, AVR and governor modules
include("SauerPaiMachineModel.jl")
include("EXST1Model.jl")
include("GasTGModel.jl")

# Use our modules
using .SauerPaiMachineModel
using .EXST1Model
using .GasTGModel

# Define constants for state indices 
# Machine state indices
const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6

# AVR state indices
const VF_IDX = 1
const VT_IDX = 2
const VLL_IDX = 3
const VFB_IDX = 4

# Governor state indices
const FV_IDX = 1
const FF_IDX = 2
const ET_IDX = 3

# Define parameters for ODE solver
mutable struct MachinePlusControllersParams
    machine::SauerPaiMachine
    avr::EXST1
    governor::GasTG
    V_terminal::Complex{Float64}  # Fixed terminal voltage
    machine_idx::UnitRange{Int64}
    avr_idx::UnitRange{Int64}
    gov_idx::UnitRange{Int64}
end

# Main function to run the simulation
function run_machine_avr_governor()
    # Define initial conditions
    V_mag = 1.0142
    V_angle = -0.007843738409384566
    P = 0.8
    Q = 0.10751144041475173

    # Convert voltage to complex form
    V_terminal = V_mag * exp(im * V_angle)
    println("Initial conditions:")
    println("V_terminal = $V_mag p.u ∠ $V_angle rad")
    println("P = $P pu, Q = $Q pu")

    # Create model instances
    machine = SauerPaiMachine(
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

    avr = EXST1(
        V_ref=V_mag  # Set reference voltage to initial voltage
    )

    governor = GasTG(
        P_ref=P  # Set reference power to initial power
    )

    # Initialize states
    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)
    avr_states = initialize_avr(avr, V_mag, Vf_init)
    ω_init = 1.0
    governor_states = initialize_gov(governor, τ_m_init, ω_init)

    # Print initial states for debugging
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

    println("Governor states:")
    println("  Fuel Value: $(governor_states[FV_IDX])")
    println("  Fuel Flow: $(governor_states[FF_IDX])")
    println("  Exhaust Temp: $(governor_states[ET_IDX])")

    println("Initial field voltage (Vf): $Vf_init")
    println("Initial mechanical torque (τm): $τ_m_init")

    # Combine states for ODE solver
    states = vcat(machine_states, avr_states, governor_states)

    # Define state indices
    machine_idx = 1:6
    avr_idx = 7:10
    gov_idx = 11:13

    # Create params for ODE solver
    p = MachinePlusControllersParams(
        machine,
        avr,
        governor,
        V_terminal,
        machine_idx,
        avr_idx,
        gov_idx
    )

    # Define auxiliary variables
    τm_aux = Float64[]
    Vf_aux = Float64[]
    ω_aux = Float64[]
    t_aux = Float64[]

    # Pack aux variables with initial values
    push!(τm_aux, τ_m_init)
    push!(Vf_aux, Vf_init)
    push!(ω_aux, ω_init)
    push!(t_aux, 0.0)

    # Machine + AVR + Governor ODE function
    function machine_avr_gov_dynamics!(du, u, params, t)
        # Extract state variables
        machine_states = u[params.machine_idx]
        avr_states = u[params.avr_idx]
        gov_states = u[params.gov_idx]

        # Prepare derivative containers
        du_machine = zeros(Float64, length(params.machine_idx))
        du_avr = zeros(Float64, length(params.avr_idx))
        du_gov = zeros(Float64, length(params.gov_idx))

        # Step 1: Update governor with current rotor speed
        τm = update_gov_states!(
            gov_states,
            du_gov,
            ω_aux[end],
            params.governor
        )

        # Step 2: Update AVR with current terminal voltage magnitude
        Vf = update_avr_states!(
            avr_states,
            du_avr,
            abs(params.V_terminal),
            params.avr
        )

        # Step 3: Update machine with current field voltage and mechanical torque
        I_terminal_machine_pos, S_terminal_machine, ω_machine, V_mag, I_mag = update_machine_states!(
            machine_states,
            du_machine,
            params.V_terminal,
            Vf_aux[end],
            τm_aux[end],
            params.machine
        )

        # Update auxiliary variables
        push!(τm_aux, τm)
        push!(Vf_aux, Vf)
        push!(ω_aux, ω_machine)
        push!(t_aux, t)

        # Copy derivatives to output
        du[params.machine_idx] .= du_machine
        du[params.avr_idx] .= du_avr
        du[params.gov_idx] .= du_gov

        # For debugging, print at integer time steps
        if abs(t - round(t)) < 0.0001
            println("t=$t: δ=$(machine_states[DELTA]), ω=$(machine_states[OMEGA]), τm=$τm, Vf=$(avr_states[VF_IDX]), V_mag=$V_mag")
        end
    end

    tspan = (0.0, 30.0)
    prob = ODEProblem(machine_avr_gov_dynamics!, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [5.0, 10.0, 15.0]

    # Define the condition for which to apply a perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end

    # Define the perturbation
    function affect!(integrator)
        #### Uncomment the desired perturbation ####
        # Load Jump
        #integrator.u[integrator.p.machine_idx[DELTA]] *= 0.85
        if isapprox(integrator.t, 5.0)
            integrator.u[integrator.p.machine_idx[DELTA]] *= 0.85
        elseif isapprox(integrator.t, 10.0)
            integrator.u[integrator.p.machine_idx[OMEGA]] *= 0.98
        elseif isapprox(integrator.t, 15.0)
            integrator.u[integrator.p.avr_idx[VF_IDX]] *= 0.90
        end

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

    delta_values = [sol[machine_idx[DELTA], i] for i in 1:length(t)]
    omega_values = [sol[machine_idx[OMEGA], i] for i in 1:length(t)]
    Vf_values = [sol[avr_idx[VF_IDX], i] for i in 1:length(t)]

    # Plot machine states
    p1 = plot(t, delta_values,
        label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, omega_values,
        label="Rotor speed (ω)", linewidth=2)

    # Plot machine fluxes
    p2 = plot(t, [sol[machine_idx[EQ_P], i] for i in 1:length(t)],
        label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[machine_idx[ED_P], i] for i in 1:length(t)],
        label="ed'", linewidth=2)
    plot!(p2, t, [sol[machine_idx[PSI_D_PP], i] for i in 1:length(t)],
        label="ψd''", linewidth=2)
    plot!(p2, t, [sol[machine_idx[PSI_Q_PP], i] for i in 1:length(t)],
        label="ψq''", linewidth=2)

    # Plot AVR states
    p3 = plot(t, Vf_values,
        label="Field Voltage", title="AVR States", linewidth=2)
    plot!(p3, t, [sol[avr_idx[VT_IDX], i] for i in 1:length(t)],
        label="Terminal Voltage", linewidth=2)
    plot!(p3, t, [sol[avr_idx[VLL_IDX], i] for i in 1:length(t)],
        label="Lead-Lag State", linewidth=2)
    plot!(p3, t, [sol[avr_idx[VFB_IDX], i] for i in 1:length(t)],
        label="Feedback State", linewidth=2)

    # Plot governor states
    p4 = plot(t, [sol[gov_idx[FV_IDX], i] for i in 1:length(t)],
        label="Fuel Value", title="Governor States", linewidth=2)
    plot!(p4, t, [sol[gov_idx[FF_IDX], i] for i in 1:length(t)],
        label="Fuel Flow", linewidth=2)
    plot!(p4, t, [sol[gov_idx[ET_IDX], i] for i in 1:length(t)],
        label="Exhaust Temp", linewidth=2)

    # Plot torques
    p5 = plot(t_aux, τm_aux,
        label="Mechanical Torque", linewidth=2)

    p_combined = plot(p1, p2, p3, p4, p5, layout=(5, 1), size=(800, 2000), left_margin=20mm)
    savefig(p_combined, "machine_avr_governor_results.png")


    return sol
end

sol = run_machine_avr_governor()

println("Simulation complete!")