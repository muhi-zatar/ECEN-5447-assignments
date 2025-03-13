# Importing the required librarires
using DifferentialEquations
using Plots
using LinearAlgebra

# Import our custom modules from other files
include("SauerPaiMachineModel.jl")
include("EXST1Model.jl")
include("GasTGModel.jl")
include("SingleMassModel.jl")
include("PowerFlowWrapper.jl")

# Use our modules
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
    machine::SauerPaiMachine
    avr::EXST1
    governor::GasTG
    shaft::SingleMass
    system_vars::SystemVars
    apply_perturbation::Bool
    perturbation_params::PerturbationParams
    machine_idx::UnitRange{Int64}
    avr_idx::UnitRange{Int64}
    gov_idx::UnitRange{Int64}
    shaft_idx::UnitRange{Int64}
end

# Main function to run the simulation
function run_simulation(network_file)
    # Step 1: Get initial conditions from power flow
    V_mag, V_angle, P, Q = get_powerflow_results(network_file)
    
    # Convert voltage to complex form
    V_terminal = V_mag * exp(im * V_angle)
    println("Initial power flow results:")
    println("V_terminal = $V_mag p.u $V_angle rad")
    println("P = $P pu, Q = $Q pu")
    
    # Step 2: Initialize models
    # Create model instances
    machine = SauerPaiMachine(
        base_power = 100.0,
        system_base_power = 100.0,
        system_base_frequency = 60.0
    )
    
    avr = EXST1(
        V_ref = V_mag  # Set reference voltage to initial voltage
    )
    
    governor = GasTG(
        P_ref = P  # Set reference power to initial power
    )
    
    shaft = SingleMass(
        H = 30.0148,
        D = 2.0,
        system_base_frequency = 60.0
    )
    
    # Initialize states
    machine_states, Vf_init = initialize_machine(machine, V_terminal, V_angle, P, Q)
    avr_states = initialize_avr(avr, V_mag, 0.0)  # Assume zero field current initially, is that right?
    governor_states = initialize_gov(governor, P)
    shaft_states = initialize_shaft(shaft, V_angle)
    
    # Combine all states into a single vector for ODE solver
    states = vcat(machine_states, avr_states, governor_states, shaft_states)
    
    # Print initial states for debugging
    println("\nInitial States:")
    println("Machine states: $machine_states")
    println("AVR states: $avr_states")
    println("Governor states: $governor_states")
    println("Shaft states: $shaft_states")
    
    # Define state indices for easier access
    machine_idx = 1:6
    avr_idx = 7:10
    gov_idx = 11:13
    shaft_idx = 14:15
    
    # Create system variables
    system_vars = SystemVars(
        V_terminal,
        V_mag,
        P,
        Q,
        1.0  # System frequency (initially at 1.0 pu)
    )
    
    # Perturbation parameters : This is an example
    perturbation_params = PerturbationParams(
        5.0,    # start_time - When perturbation begins
        5.1,    # end_time - When perturbation ends
        0.8,    # magnitude_factor - For voltage perturbation
        1.2,    # load_factor - For load perturbation
        0.7     # fault_factor - For fault perturbation (less severe for stability)
    )
    
    # Create ODE parameters structure
    p = ODEParams(
        machine,
        avr,
        governor,
        shaft,
        system_vars,
        false,  # apply_perturbation
        perturbation_params,
        machine_idx,
        avr_idx,
        gov_idx,
        shaft_idx
    )
    
    # Define global variables for intermediate values
    # These are needed because we can't modify p during integration
    τe_global = Ref(0.0)
    τm_global = Ref(P)
    Vf_global = Ref(Vf_init)
    delta_global = Ref(V_angle)
    omega_global = Ref(1.0)
    
    # Function to apply perturbation
    function apply_perturbation!(vars, t, params)
        if t >= params.start_time && t <= params.end_time
            vars.V_terminal = vars.V_terminal * params.fault_factor
            vars.V_terminal_magnitude = abs(vars.V_terminal)
            return true
        end
        return false
    end
    
    # ODE function
    function ode_system!(du, u, p, t)
        # Extract states for each component
        machine_states = u[p.machine_idx]
        avr_states = u[p.avr_idx]
        gov_states = u[p.gov_idx]
        shaft_states = u[p.shaft_idx]
        
        # Extract parameters
        machine = p.machine
        avr = p.avr
        governor = p.governor
        shaft = p.shaft
        system_vars = p.system_vars
        
        # Apply perturbation if specified
        if p.apply_perturbation
            apply_perturbation!(system_vars, t, p.perturbation_params)
        end
        
        # Make copies of states for Float64 compatibility if needed
        machine_states_f64 = convert.(Float64, machine_states)
        avr_states_f64 = convert.(Float64, avr_states)
        gov_states_f64 = convert.(Float64, gov_states)
        shaft_states_f64 = convert.(Float64, shaft_states)
        
        # Arrays for derivatives
        du_machine = zeros(Float64, length(p.machine_idx))
        du_avr = zeros(Float64, length(p.avr_idx))
        du_gov = zeros(Float64, length(p.gov_idx))
        du_shaft = zeros(Float64, length(p.shaft_idx))
        
        # 1. Update shaft states
        delta, omega = update_shaft_states!(
            shaft_states_f64,
            du_shaft,
            τe_global[],
            τm_global[],
            system_vars.ω_sys,
            shaft
        )
        
        # Update global values
        delta_global[] = delta
        omega_global[] = omega
        
        # 2. Update machine states
        τe, I_grid = update_machine_states!(
            machine_states_f64,
            du_machine,
            delta,
            omega,
            Vf_global[],
            system_vars.V_terminal,
            machine
        )
        
        # Update global torque
        τe_global[] = τe
        
        # Update field current (approximation)
        Ifd = abs(machine_states_f64[3]) / (machine.Xd - machine.Xd_p)
        
        # 3. Update AVR states
        Vf = update_avr_states!(
            avr_states_f64,
            du_avr,
            system_vars.V_terminal_magnitude,
            0.0,  # No PSS output
            Ifd,
            avr
        )
        
        # Update global field voltage
        Vf_global[] = Vf
        
        # 4. Update governor states
        τm = update_gov_states!(
            gov_states_f64,
            du_gov,
            omega,
            governor
        )
        
        # Update global mechanical torque
        τm_global[] = τm
        
        # Copy derivatives to output
        du[p.machine_idx] .= du_machine
        du[p.avr_idx] .= du_avr
        du[p.gov_idx] .= du_gov
        du[p.shaft_idx] .= du_shaft
        
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
    combined_plot = plot(p1, p2, p3, p4, p5, p6, layout=(6,1), size=(800, 1200))
    display(combined_plot)
    savefig(combined_plot, "haha.png")
    
    return sol
end

# Run the simulation
sol = run_simulation("../data/ThreeBusNetwork.raw")  # Adjust the path to your network file

println("Simulation complete!")