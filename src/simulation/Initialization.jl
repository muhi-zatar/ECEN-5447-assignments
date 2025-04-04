module Initialization

export initialize_system

using ..Types
using ..MachineComponents
using ..AVRComponents
using ..GovernorComponents
using ..NetworkComponents
using ..PowerFlowComponents
using ..ComponentIndex
using ..Logging

"""
    initialize_system(model::PowerSystemModel, network_file::String, logger=nothing)

Initialize all component states based on power flow solution.
Updates the model with initialized states and returns the global state vector.
"""
function initialize_system(model::PowerSystemModel, network_file::String, logger=nothing)
    if logger !== nothing
        log_info(logger, "Initializing system from $network_file")
    end
    
    # Create state indices
    indices = create_state_indices(model)
    
    # Create initial global state vector
    states = zeros(Float64, indices.num_states)
    
    # Step 1: Get initial conditions from power flow
    V_sol, θ_sol, P_sol, Q_sol = get_powerflow_results(network_file)
    
    # Extract bus 2 (generator bus) values
    V_mag = V_sol[Int(BUS_MACHINE_MODEL)]
    V_angle = θ_sol[Int(BUS_MACHINE_MODEL)]
    P = P_sol[Int(BUS_MACHINE_MODEL)]
    Q = Q_sol[Int(BUS_MACHINE_MODEL)]
    
    # Convert voltage to complex form
    V_terminal_init = V_mag * exp(im * V_angle)
    I_terminal_init = conj(complex(P, Q) / V_terminal_init)
    
    if logger !== nothing
        log_info(logger, "Power flow results:")
        log_info(logger, "V_terminal = $V_mag p.u $V_angle rad")
        log_info(logger, "P = $P pu, Q = $Q pu")
        log_info(logger, "I_terminal = $(abs(I_terminal_init)) p.u., $(angle(I_terminal_init)) rad")
    end
    
    # Step 2: Initialize network
    network_states, i_2_d, i_2_q = initialize_network_states(model.network, V_sol, θ_sol, P_sol, Q_sol)
    
    # Step 3: Initialize machine if enabled
    if model.enabled_components["machine"]
        machine_states, Vf_init, τ_m_init = initialize_machine_states(
            model.machine, V_terminal_init, V_angle, P, Q
        )
        println("Machine states initialized.")
        println(machine_states)
        # Step 4: Initialize AVR if enabled
        if model.enabled_components["avr"]
            avr_states = initialize_avr_states(model.avr, V_mag, Vf_init)
        else
            # Create constant field voltage states
            avr_states = zeros(Float64, 4)
            avr_states[EFD_IDX] = Vf_init
        end
        
        # Step 5: Initialize governor if enabled
        if model.enabled_components["governor"]
            governor_states = initialize_gov_states(model.governor, τ_m_init, 1.0)
        else
            # Create constant mechanical torque states
            governor_states = zeros(Float64, 3)
            governor_states[FV_IDX] = τ_m_init
            governor_states[FF_IDX] = τ_m_init
            governor_states[ET_IDX] = τ_m_init
        end
    else
        # If machine is disabled, use default values
        machine_states = zeros(Float64, 8)
        machine_states[OMEGA] = 1.0  # Set speed to nominal
        
        # Simple AVR states
        avr_states = zeros(Float64, 4)
        avr_states[EFD_IDX] = 1.0  # Default field voltage
        
        # Simple governor states
        governor_states = zeros(Float64, 3)
        governor_states[FV_IDX] = 1.0  # Default mechanical torque
        governor_states[FF_IDX] = 1.0
        governor_states[ET_IDX] = 1.0
    end
    
    # Combine all states into the global state vector
    states[indices.network_range] = network_states
    states[indices.machine_range] = machine_states
    states[indices.avr_range] = avr_states
    states[indices.governor_range] = governor_states
    
    if logger !== nothing
        log_info(logger, "System initialized with $(length(states)) states")
    end
    
    return states
end

end # module