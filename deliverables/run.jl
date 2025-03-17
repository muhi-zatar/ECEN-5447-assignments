# Ensure script runs in the correct environment
cd(@__DIR__)
using Pkg
Pkg.activate("../.")
#Pkg.resolve()
#Pkg.instantiate()

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
    network_states, i_2_r_init, i_2_i_init = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)
    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)
    avr_states = initialize_avr(avr, V_mag, Vf_init)  # Assume zero field current initially, is that right?
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
    Vf_global = Vf_init
    τm_global = P
    V_terminal_magnitude_global = V_mag
    V_terminal_global = V_terminal
    θ_global = V_angle
    ω_global = 1.0
    i_2_r_global = i_2_r_init                    # Current injection at Bus 2 (network D-axis)
    i_2_i_global = i_2_i_init                    # Current injection at Bus 2 (network Q-axis)

    # Function to apply perturbation
    # function apply_perturbation!(vars, t, params)
    #     if t >= params.start_time && t <= params.end_time
    #         vars.V_terminal = vars.V_terminal * params.fault_factor
    #         vars.V_terminal_magnitude = abs(vars.V_terminal)
    #         return true
    #     end
    #     return false
    # end

    # Step 2.5: Define Mass Matrix for the system
    M_system = zeros(Float64, 35)
    M_system[network_idx] = network.M
    M_system[machine_idx] = ones(6)
    M_system[avr_idx] = ones(4)
    M_system[gov_idx] = ones(3)

    # ODE function
    function ode_system!(du, u, p, t, mass_matrix=M_system)
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
        Vf_avr = update_avr_states!(avr_states_f64, du_avr, V_terminal_magnitude_global, avr)

        # Update global field voltage
        Vf_global = Vf_avr

        τ_m = update_gov_states!(gov_states_f64, du_gov, ω_global, governor)

        # Update global mechanical torque
        τm_global = τ_m

        V_mag_bus, I_RI, θ_bus, ω_machine = update_machine_states!(machine_states_f64, du_machine, V_terminal_global, Vf_global, τm_global, machine)

        θ_global = θ_bus

        i_2_r_global = real(I_RI)
        i_2_i_global = imag(I_RI)

        V_terminal_network = update_network_states!(network_states_f64, du_network, i_2_r_global, i_2_i_global, θ_global, network)

        # Update global terminal voltage and rotor speed
        ω_global = ω_machine
        V_terminal_magnitude_global = V_mag_bus             # Bus 2 terminal voltage magnitude
        # TODO: Need to make the proper conversion from complex terminal voltage in the network DQ 
        # frame to the positive sequence domain.
        # Alternatively, change update_machine_states! to only use V_dq (I'm leaning towards this)
        V_terminal_global = V_terminal_network                             # Bus 2 terminal voltage (positive sequence phasor)


        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.machine_idx] .= du_machine
        du[p.avr_idx] .= du_avr
        du[p.gov_idx] .= du_gov

        # For debugging printing some results,
        if abs(t - round(t)) < 0.0001
            #println("t=$t: delta=$delta, omega=$omega, τe=$τe, τm=$τm, Vf=$Vf, V_mag=$(system_vars.V_terminal_magnitude)")
            println("t=$t: omega=$ω_global, τm=$τm_global, Vf=$Vf_global, V_mag=$V_terminal_magnitude_global")
        end
    end

    # Step 4: Set up and solve the ODE system
    tspan = (0.0, 5.0)
    prob = ODEProblem(ode_system!, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [4.5]               # Setting this far ahead for now – we can change this later

    # Define the condition for which to apply a perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end

    # Define the perturbation
    function affect!(integrator)
        #### Uncomment the desired perturbation ####
        # Load Jump
        integrator.p.network.Z_L *= 1.15

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

    # Process results
    t = sol.t
    #voltage_magnitudes = fill(NaN, length(t))

    # for (i, ti) in enumerate(t)
    #     # Calculate voltage at each saved time point
    #     if ti >= perturbation_params.start_time && ti <= perturbation_params.end_time
    #         voltage_magnitudes[i] = V_mag * perturbation_params.fault_factor
    #     else
    #         voltage_magnitudes[i] = V_mag
    #     end
    # end

    # Create plots - Just as an example, 90% of them are not needed
    p1 = plot(t, [sol[machine_idx[3], i] for i in 1:length(t)], label="eq'", title="Machine Fluxes")
    plot!(p1, t, [sol[machine_idx[4], i] for i in 1:length(t)], label="ed'")
    savefig(p1, "fluxes.png")
    p2 = plot(t, [sol[avr_idx[1], i] for i in 1:length(t)], label="Field Voltage", title="AVR")
    savefig(p2, "vf.png")
    #p3 = plot(t, [sol[gov_idx[2], i] for i in 1:length(t)], label="Mechanical Power", title="Governor")
    p3 = plot(t, [sol[network_idx[9], i] for i in 1:length(t)], title="Terminal Voltage (Machine DQ)", label="Vd", linewidth=2)
    plot!(p3, t, [sol[network_idx[10], i] for i in 1:length(t)], label="Vq", linewidth=2)
    savefig(p3, "vt.png")

    # Combine plots
    #combined_plot = plot(p1, p2, p3, layout=(3, 1), size=(800, 1200))
    #display(combined_plot)
    #savefig(combined_plot, "haha.png")

    return sol
end

# Run the simulation
sol = run_simulation("../data/ThreeBusMultiLoad.raw")  # Adjust the path to your network file

println("Simulation complete!")