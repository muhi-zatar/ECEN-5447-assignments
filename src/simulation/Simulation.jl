module Simulation

export ODEParams, run_simulation

using DifferentialEquations
using LinearAlgebra
using ..Types
using ..ComponentIndex: create_state_indices, StateIndices
using ..MachineComponents
using ..AVRComponents
using ..GovernorComponents
using ..NetworkComponents
using ..Initialization
using ..Perturbation: apply_perturbation!, PerturbationType, LINE_TRIP, LOAD_INCREASE, LOAD_DECREASE
using ..Logging
using ..Plotting
using ..Transformations: ri_dq, dq_ri

"""
    ODEParams

Contains all parameters needed for the ODE solver.
"""
struct ODEParams
    model::PowerSystemModel
    indices::StateIndices
    logger::Any
end

"""
    run_simulation(model::PowerSystemModel, logger=nothing)

Run the power system dynamics simulation.
Returns the simulation solution.
"""
function run_simulation(model::PowerSystemModel, logger=nothing)
    # Create indices
    indices = create_state_indices(model)
    
    # Initialize components
    states = zeros(Float64, indices.num_states)
    network_path = "../../data/ThreeBusMultiLoad.raw"
    states = initialize_system(model, network_path, logger)

    # Extract simulation parameters
    tspan = model.simulation_params.tspan
    dt = model.simulation_params.dt
    save_interval = model.simulation_params.save_interval
    perturb_time = model.simulation_params.perturb_time
    
    # Create ODE parameters structure
    p = ODEParams(model, indices, logger)
    
    # Build mass matrix (diagonal) for DAE system
    M_system = zeros(Float64, indices.num_states)
    
    # Set network mass matrix components
    if isa(model.network.M, Vector)
        M_system[indices.network_range] .= model.network.M
    else
        for i in 1:length(indices.network_range)
            M_system[indices.network_range[i]] = model.network.M[i, i]
        end
    end
    
    # Set machine, AVR, and governor mass matrix components
    M_system[indices.machine_range] .= 1.0
    M_system[indices.avr_range] .= 1.0
    M_system[indices.governor_range] .= 1.0
    
    # Create mass matrix for the DAE solver
    mass_matrix = Diagonal(M_system)

    # Define ODE system function
    function ode_system!(du, u, p, t)
        # Extract parameters
        model = p.model
        indices = p.indices
        logger = p.logger
        
        # Extract component states
        network_states = view(u, indices.network_range)
        machine_states = view(u, indices.machine_range)
        avr_states = view(u, indices.avr_range)
        gov_states = view(u, indices.governor_range)
        
        # Extract derivatives for each component
        du_network = view(du, indices.network_range)
        du_machine = view(du, indices.machine_range)
        du_avr = view(du, indices.avr_range)
        du_gov = view(du, indices.governor_range)

        # Get terminal voltage from network
        v_2_d = network_states[indices.network[:V_2_D_IDX] - indices.network_range.start + 1]
        v_2_q = network_states[indices.network[:V_2_Q_IDX] - indices.network_range.start + 1]
        V_RI = dq_ri(0.0) * [v_2_d; v_2_q]  # Always use 0.0 for network reference frame
        V_terminal = complex(V_RI[1], V_RI[2])
        V_terminal_mag = abs(V_terminal)
        
        # Get machine speed
        ω = machine_states[indices.machine[:OMEGA] - indices.machine_range.start + 1]
        
        # Get field voltage from AVR
        efd = avr_states[indices.avr[:EFD_IDX] - indices.avr_range.start + 1]

        # Update AVR states if enabled
        if model.enabled_components["avr"]
            update_avr_states!(avr_states, du_avr, V_terminal_mag, model.avr)
        else
            # Use constant field voltage
            du_avr .= 0.0
        end
        
        # Update governor states if enabled
        if model.enabled_components["governor"]
            τ_m = update_gov_states!(gov_states, du_gov, ω, model.governor)
        else
            # Use constant mechanical torque (from initialization)
            τ_m = gov_states[indices.governor[:FV_IDX] - indices.governor_range.start + 1]
            du_gov .= 0.0
        end
        # # Update machine states if enabled
        if model.enabled_components["machine"]
            _, S_terminal_machine, _, _, _ = update_machine_states!(
                machine_states, du_machine, V_terminal, efd, τ_m, model.machine
            )
        else
            du_machine .= 0.0
            # Calculate power from terminal voltage and machine current
            i_d = 0.0  # Would need to define a constant current source model
            i_q = 0.0
            S_terminal_machine = complex(v_2_d * i_d + v_2_q * i_q, v_2_q * i_d - v_2_d * i_q)
        end
        
        # Update network states
        _, _, _ = update_network_states!(network_states, du_network, S_terminal_machine, model.network)
        
        # Log state values if logger is provided
        if logger !== nothing && mod(t, 0.1) < dt
            log_states(logger, t, u, indices)
        end
    end
    
    # Create ODE function with mass matrix
    dae_problem = ODEFunction(ode_system!, mass_matrix=mass_matrix)
    
    # Define callback for perturbation
    perturb_times = [perturb_time]
    
    # Define condition for perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end
    
    # Define perturbation effect
    function affect!(integrator)
        # Apply perturbation to the model
        apply_perturbation!(
            integrator.p.model,
            perturbation_type=LOAD_DECREASE,
            parameter_changes=Dict{Symbol, Any}()
        )
                
        # Log the perturbation
        if integrator.p.logger !== nothing
            log_event(integrator.p.logger, integrator.t, "Perturbation applied")
        end
    end
    
    # Create callback
    cb = DiscreteCallback(condition, affect!)
    
    # Set up ODE problem
    prob = ODEProblem(dae_problem, states, tspan, p)
    
    # Solve the ODE system
    if logger !== nothing
        log_info(logger, "Starting simulation with $(indices.num_states) states")
    end
    
    # Select appropriate solver for DAE systems
    sol = solve(
        prob, 
        Rodas5(autodiff=false), 
        dt=dt, 
        adaptive=false, 
        saveat=save_interval, 
        callback=cb, 
        tstops=perturb_times
    )
    
    # Log simulation complete
    if logger !== nothing
        log_info(logger, "Simulation complete: $(length(sol.t)) time points")
    end
    
    return sol
end

end # module