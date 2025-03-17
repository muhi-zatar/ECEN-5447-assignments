# Ensure script runs in the correct environment
cd(@__DIR__)
using Pkg
Pkg.activate("../.")
Pkg.resolve()
Pkg.instantiate()

# Importing the required librarires
using DifferentialEquations
using Plots
using LinearAlgebra

# Import our custom modules from other files
include("NetworkModel.jl")
include("SauerPaiMachineModel.jl")
include("EXST1Model.jl")
include("GasTGModel.jl")
include("SingleMassModel.jl")
include("PowerFlowWrapper.jl")

# Use our modules
using .NetworkModel
using .SauerPaiMachineModel
using .EXST1Model
using .GasTGModel
using .SingleMassModel
using .PowerFlowWrapper

# System variables structure
mutable struct SystemVars
    V_terminal::Complex{Float64}
    V_terminal_magnitude::Float64
    P_load::Float64
    Q_load::Float64
    ω_sys::Float64
end

# Perturbation parameters
mutable struct PerturbationParams
    start_time::Float64
    end_time::Float64
    magnitude_factor::Float64
    load_factor::Float64
    fault_factor::Float64
end

# Define parameters for ODE solver
struct ODEParams
    network::ThreeBusNetwork
    machine::SauerPaiMachine
    avr::EXST1
    governor::GasTG
    network_idx::UnitRange{Int64}
    machine_idx::UnitRange{Int64}
    avr_idx::UnitRange{Int64}
    gov_idx::UnitRange{Int64}
end

# Main function to run the simulation
function run_simulation(network_file)
    # Step 1: Get initial conditions from power flow (all buses)
    V_sol, θ_sol, P_sol, Q_sol = get_powerflow_results(network_file)

    # Select bus for our model :D
    V_mag = V_sol[BUS_MACHINE_MODEL]
    V_angle = θ_sol[BUS_MACHINE_MODEL]
    P = P_sol[BUS_MACHINE_MODEL]
    Q = Q_sol[BUS_MACHINE_MODEL]

    # Convert voltage to complex form
    V_terminal = V_mag * exp(im * V_angle)
    println("Initial power flow results:")
    println("V_terminal = $V_mag p.u $V_angle rad")
    println("P = $P pu, Q = $Q pu")

    # Step 2: Initialize models
    # Create model instances
    network = ThreeBusNetwork()

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
    network_states, i_2_d_init, i_2_q_init = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)
    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)
    avr_states, Vs = initialize_avr(avr, V_mag, Vf_init)  # Assume zero field current initially, is that right?
    ω_init = 1.0
    governor_states = initialize_gov(governor, τ_m_init, ω_init)

    # Combine all states into a single vector for ODE solver
    states = vcat(network_states, machine_states, avr_states, governor_states)

    # Print initial states for debugging
    println("\nInitial States:")
    println("Machine states: $machine_states")
    println("AVR states: $avr_states")
    println("Governor states: $governor_states")

    # Define state indices for easier access
    network_idx = 1:22
    machine_idx = 23:28
    avr_idx = 29:32
    gov_idx = 33:35

    # Create ODE parameters structure
    p = ODEParams(
        network,
        machine,
        avr,
        governor,
        network_idx,
        machine_idx,
        avr_idx,
        gov_idx,
    )

    # Define global variables for intermediate values
    # These are needed because we can't modify p during integration
    Vf_global = Ref(Vf_init)
    τm_global = Ref(P)
    V_terminal_magnitude_global = Ref(V_mag)
    V_terminal_global = Ref(V_terminal)
    ω_global = Ref(1.0)
    i_2_d_global = Ref(i_2_d_init)                    # Current injection at Bus 2 (network D-axis)
    i_2_q_global = Ref(i_2_q_init)                    # Current injection at Bus 2 (network Q-axis)

    # Function to apply perturbation
    # function apply_perturbation!(vars, t, params)
    #     if t >= params.start_time && t <= params.end_time
    #         vars.V_terminal = vars.V_terminal * params.fault_factor
    #         vars.V_terminal_magnitude = abs(vars.V_terminal)
    #         return true
    #     end
    #     return false
    # end

    # ODE function
    function ode_system!(du, u, p, t)
        # Extract states for each component
        network_states = u[p.network_idx]
        machine_states = u[p.machine_idx]
        avr_states = u[p.avr_idx]
        gov_states = u[p.gov_idx]

        # Extract parameters
        network = p.network
        machine = p.machine
        avr = p.avr
        governor = p.governor

        # Apply perturbation if specified
        # if p.apply_perturbation
        #     apply_perturbation!(system_vars, t, p.perturbation_params)
        # end

        # Make copies of states for Float64 compatibility if needed
        network_states_f64 = convert.(Float64, network_states)
        machine_states_f64 = convert.(Float64, machine_states)
        avr_states_f64 = convert.(Float64, avr_states)
        gov_states_f64 = convert.(Float64, gov_states)

        # Arrays for derivatives
        du_network = zeros(Float64, length(p.network_idx))
        du_machine = zeros(Float64, length(p.machine_idx))
        du_avr = zeros(Float64, length(p.avr_idx))
        du_gov = zeros(Float64, length(p.gov_idx))

        # Update the states of each component
        update_avr_states!(avr_states_f64, du_avr, V_terminal_magnitude_global, avr)

        # Update global field voltage
        Vf_global[] = avr_states_f64[1]

        τ_m = update_gov_states!(gov_states_f64, du_gov, ω_global, governor)

        # Update global mechanical torque
        τm_global[] = τ_m

        V_dq, I_dq = update_machine_states!(machine_states_f64, du_machine, V_terminal_global, Vf_global, τm_global, machine)

        # Update global bus 2 current injection
        i_2_d_global[] = I_dq[1]
        i_2_q_global[] = I_dq[2]

        update_network_states!(network_states_f64, du_network, i_2_d_global, i_2_q_global, network)

        # Update global terminal voltage and rotor speed
        ω_global[] = ω
        V_terminal_magnitude_global[] = abs(V_dq)             # Bus 2 terminal voltage magnitude
        # TODO: Need to make the proper conversion from complex terminal voltage in the network DQ 
        # frame to the positive sequence domain.
        # Alternatively, change update_machine_states! to only use V_dq (I'm leaning towards this)
        V_terminal_global[] = 0.0                             # Bus 2 terminal voltage (positive sequence phasor)


        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.machine_idx] .= du_machine
        du[p.avr_idx] .= du_avr
        du[p.gov_idx] .= du_gov

        # For debugging printing some results,
        if abs(t - round(t)) < 0.0001
            println("t=$t: delta=$delta, omega=$omega, τe=$τe, τm=$τm, Vf=$Vf, V_mag=$(system_vars.V_terminal_magnitude)")
        end
    end

    # Step 4: Set up and solve the ODE system
    tspan = (0.0, 20.0)
    prob = ODEProblem(ode_system!, states, tspan, p)

    # Using Euler or other simple explicit method that doesn't use automatic differentiation
    # Facing an error with autmatic differentiation
    # TODO: check why
    sol = solve(prob, BS3(), dt=0.00005, adaptive=false, saveat=0.01)

    # Process results
    t = sol.t
    voltage_magnitudes = fill(NaN, length(t))

    for (i, ti) in enumerate(t)
        # Calculate voltage at each saved time point
        if ti >= perturbation_params.start_time && ti <= perturbation_params.end_time
            voltage_magnitudes[i] = V_mag * perturbation_params.fault_factor
        else
            voltage_magnitudes[i] = V_mag
        end
    end

    # Create plots - Just as an example, 90% of them are not needed
    p1 = plot(t, [sol[shaft_idx[1], i] for i in 1:length(t)], label="Rotor Angle (rad)", title="Rotor Dynamics")
    p2 = plot(t, [sol[shaft_idx[2], i] for i in 1:length(t)], label="Rotor Speed (pu)", title="Speed")
    p3 = plot(t, [sol[machine_idx[3], i] for i in 1:length(t)], label="eq'", title="Machine Fluxes")
    plot!(p3, t, [sol[machine_idx[4], i] for i in 1:length(t)], label="ed'")
    p4 = plot(t, [sol[avr_idx[3], i] for i in 1:length(t)], label="Regulator Output", title="AVR")
    p5 = plot(t, [sol[gov_idx[2], i] for i in 1:length(t)], label="Mechanical Power", title="Governor")
    p6 = plot(t, voltage_magnitudes, label="Terminal Voltage (pu)", title="Voltage", linewidth=2)

    # Combine plots
    combined_plot = plot(p1, p2, p3, p4, p5, p6, layout=(6, 1), size=(800, 1200))
    display(combined_plot)
    savefig(combined_plot, "haha.png")

    return sol
end

# Run the simulation
sol = run_simulation("../data/ThreeBusNetwork.raw")  # Adjust the path to your network file

println("Simulation complete!")