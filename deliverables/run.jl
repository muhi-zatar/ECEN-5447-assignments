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
    I_terminal = conj(complex(P, Q) / V_terminal)
    println("Initial power flow results:")
    println("V_terminal = $V_mag p.u $V_angle rad")
    println("P = $P pu, Q = $Q pu")
    println("I_terminal = $(abs(I_terminal)) p.u., $(angle(I_terminal)) rad")

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
    network_idx = 1:16
    machine_idx = 17:22
    avr_idx = 23:26
    gov_idx = 27:29

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

    # Define auxiliary variables for intermediate values.
    # We will use vectors rather than global variables to avoid performance issues (https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-global-variables)
    # These are needed because we can't modify p during integration
    # Start by allocating empty typed vectors
    # Used by ODE System:
    V_terminal_magnitude_aux = Float64[]          # Machine bus voltage magnitude (from Power Flow)
    ω_aux = Float64[]                             # Rotor angular velocity (assumed to be 1.0 p.u.)
    V_terminal_aux = Complex{Float64}[]           # Machine bus quasi-static voltage phasor (from Power Flow)
    Vf_aux = Float64[]                            # Field voltage (from machine initialization)
    τm_aux = Float64[]                            # Mechanical torque (Assume τm = τe = P)
    S_terminal_machine_aux = Complex{Float64}[]   # Machine bus apparent power (from Power Flow)

    # For debugging:
    S_terminal_network_aux = Complex{Float64}[]
    V_mag_machine_aux = Float64[]
    I_mag_machine_aux = Float64[]
    t_aux = Float64[]

    # Now populate the first entry with the appropriate values from initialization
    push!(V_terminal_magnitude_aux, V_mag)
    push!(ω_aux, ω_init)
    push!(V_terminal_aux, V_terminal)
    push!(Vf_aux, Vf_init)
    push!(τm_aux, τ_m_init)
    push!(S_terminal_machine_aux, Complex(P, Q))

    # For debugging:
    push!(S_terminal_network_aux, Complex(0.0, 0.0))
    push!(V_mag_machine_aux, 0.0)
    push!(I_mag_machine_aux, 0.0)
    push!(t_aux, 0.0)
    println("Initial auxiliary variables:\nomega=$(ω_aux[end]), τm=$(τm_aux[end]), Vf=$(Vf_aux[end]), V_mag=$(V_terminal_magnitude_aux[end]), V_mag_machine=$(V_mag_machine_aux[end]), I_mag_machine=$(I_mag_machine_aux[end]), S_terminal_machine=$(S_terminal_machine_aux[end]), S_terminal_network=$(S_terminal_network_aux[end])")

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
    # M_system = zeros(Float64, 35)
    # M_system[network_idx] = network.M
    # M_system[machine_idx] = ones(6)
    # M_system[avr_idx] = ones(4)
    # M_system[gov_idx] = ones(3)

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
        efd = update_avr_states!(avr_states_f64, du_avr, V_terminal_magnitude_aux[end], avr)

        τ_m = update_gov_states!(gov_states_f64, du_gov, ω_aux[end], governor)

        I_terminal_machine_pos, S_terminal_machine, ω_machine, V_mag, I_mag = update_machine_states!(machine_states_f64, du_machine, V_terminal_aux[end], Vf_aux[end], τm_aux[end], machine)

        V_terminal_network_pos, S_terminal_network, I_terminal_network_pos, i_2_d, i_2_q = update_network_states!(network_states_f64, du_network, S_terminal_machine_aux[end], network)

        # Update global/auxiliary variables to prepare for the next step
        # Machine bus voltage magnitude: Used by AVR, returned by network
        push!(V_terminal_magnitude_aux, abs(V_terminal_network_pos))

        # Rotor angular velocity: Used by governor, returned by machine
        push!(ω_aux, ω_machine)

        # Machine bus quasi-static voltage phasor: Used by machine, returned by network
        push!(V_terminal_aux, V_terminal_network_pos)

        # Field voltage: Used by machine, returned by AVR
        push!(Vf_aux, efd)

        # Mechanical torque: Used by machine, returned by governor
        push!(τm_aux, τ_m)

        # Machine apparent power: Used by network, returned by machine
        push!(S_terminal_machine_aux, S_terminal_machine)

        # Debugging:
        push!(S_terminal_network_aux, S_terminal_network)
        push!(V_mag_machine_aux, V_mag)
        push!(I_mag_machine_aux, I_mag)
        push!(t_aux, t)


        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.machine_idx] .= du_machine
        du[p.avr_idx] .= du_avr
        du[p.gov_idx] .= du_gov

        # For debugging printing some results,
        if abs(t - round(t)) < 0.000001
            println("t=$t")
            #println("t=$t: omega=$(ω_aux[end]), τm=$(τm_aux[end]), Vf=$(Vf_aux[end]), V_mag=$(V_terminal_magnitude_aux[end]), V_mag_machine=$(V_mag_machine_aux[end]), I_mag_machine=$(I_mag_machine_aux[end]), S_terminal_machine=$(S_terminal_machine_aux[end]), S_terminal_network=$(S_terminal_network_aux[end])")
            println("Machine States: δ=$(machine_states_f64[1]), ω=$(machine_states_f64[2]), EQ_P=$(machine_states_f64[3]), ED_P=$(machine_states_f64[4]), PSI_D_PP=$(machine_states_f64[5]), PSI_Q_PP=$(machine_states_f64[6])")
            println("AVR States: VF=$(avr_states_f64[1]), VT=$(avr_states_f64[2]), VLL=$(avr_states_f64[3]), VFB=$(avr_states_f64[4])")
            println("Gov States: FV=$(gov_states_f64[1]), FF=$(gov_states_f64[2]), ET=$(gov_states_f64[3])")
            #println("Network States: I_12_D=$(network_states_f64[1]), I_12_Q=$(network_states_f64[2]), I_13_D=$(network_states_f64[3]), I_13_Q=$(network_states_f64[4]), V_2_D=$(network_states_f64[9]), V_2_Q=$(network_states_f64[10]), I_B2_D=$(network_states_f64[19]), I_B2_Q=$(network_states_f64[20])")
        end
    end

    # Step 4: Set up and solve the ODE system
    tspan = (0.0, 10.0)
    prob = ODEProblem(ode_system!, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [15.0]               # Setting this far ahead for now – we can change this later

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
    delta_values = [sol[machine_idx[1], i] for i in 1:length(t)]
    omega_values = [sol[machine_idx[2], i] for i in 1:length(t)]
    Vf_values = [sol[avr_idx[1], i] for i in 1:length(t)]

    # Create plots - Just as an example, 90% of them are not needed
    # Plot shaft states
    p1 = plot(t, delta_values,
        label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, omega_values,
        label="Rotor speed (ω)", linewidth=2)
    savefig(p1, "shaft_states.png")

    # Plot fluxes
    p2 = plot(t, [sol[machine_idx[3], i] for i in 1:length(t)],
        label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[machine_idx[4], i] for i in 1:length(t)],
        label="ed'", linewidth=2)
    plot!(p2, t, [sol[machine_idx[5], i] for i in 1:length(t)],
        label="ψd''", linewidth=2)
    plot!(p2, t, [sol[machine_idx[6], i] for i in 1:length(t)],
        label="ψq''", linewidth=2)
    savefig(p2, "fluxes.png")

    # Plot avr states
    p3 = plot(t, Vf_values,
        label="Field Voltage", title="AVR States", linewidth=2)
    plot!(p3, t, [sol[avr_idx[2], i] for i in 1:length(t)],
        label="Terminal Voltage", linewidth=2)
    plot!(p3, t, [sol[avr_idx[3], i] for i in 1:length(t)],
        label="Lead-Lag State", linewidth=2)
    plot!(p3, t, [sol[avr_idx[4], i] for i in 1:length(t)],
        label="Feedback State", linewidth=2)
    savefig(p3, "avr_states.png")

    # Plot governor states
    p4 = plot(t, [sol[gov_idx[1], i] for i in 1:length(t)],
        label="Fuel Value", title="Governor States", linewidth=2)
    plot!(p4, t, [sol[gov_idx[2], i] for i in 1:length(t)],
        label="Fuel Flow", linewidth=2)
    plot!(p4, t, [sol[gov_idx[3], i] for i in 1:length(t)],
        label="Exhaust Temp", linewidth=2)
    savefig(p4, "gov_states.png")

    # Plot torques
    p5 = plot(t_aux, τm_aux,
        label="Mechanical Torque", title="Machine Torques", linewidth=2)
    savefig(p5, "torque.png")

    # Plot voltages
    p6 = plot(t_aux, V_terminal_magnitude_aux, label="Terminal Voltage Magnitude", title="Voltages", linewidth=2)
    plot!(p6, t_aux, Vf_aux, label="Field Voltage", linewidth=2)
    savefig(p6, "voltages.png")

    # Plot powers
    p7 = plot(t_aux, real(S_terminal_machine_aux), label="P", title="Power Injections at Machine Bus", linewidth=2)
    plot!(p7, t_aux, imag(S_terminal_machine_aux), label="Q", linewidth=2)
    savefig(p7, "power.png")

    # p_combined = plot(p1, p2, p3, p4, p5, p6, p7, layout=(7, 1), size=(1000, 2400))
    # savefig(p_combined, "machine_avr_governor_results.png")

    return sol
end

# Run the simulation
sol = run_simulation("../data/ThreeBusMultiLoad.raw")  # Adjust the path to your network file

println("Simulation complete!")