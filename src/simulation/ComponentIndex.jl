module ComponentIndex

export StateIndices, create_state_indices

using ..MachineComponents
using ..AVRComponents
using ..GovernorComponents
using ..NetworkComponents

"""
    StateIndices

Stores the indices for all components in the global state vector.
This allows components to know where their states are located.
"""
struct StateIndices
    # Network state indices with their global positions
    network::Dict{Symbol, Int}
    
    # Machine state indices with their global positions
    machine::Dict{Symbol, Int}
    
    # AVR state indices with their global positions
    avr::Dict{Symbol, Int}
    
    # Governor state indices with their global positions
    governor::Dict{Symbol, Int}
    
    # Component ranges in the global state vector
    network_range::UnitRange{Int64}
    machine_range::UnitRange{Int64}
    avr_range::UnitRange{Int64}
    governor_range::UnitRange{Int64}
    
    # Total number of states
    num_states::Int
end

"""
    create_state_indices(model)

Creates a StateIndices object for a given power system model.
"""
function create_state_indices(model)
    # Network constants
    network_const = [
        :I_12_D_IDX, :I_12_Q_IDX, 
        :I_13_D_IDX, :I_13_Q_IDX, 
        :I_23_D_IDX, :I_23_Q_IDX, 
        :V_1_D_IDX, :V_1_Q_IDX, 
        :V_2_D_IDX, :V_2_Q_IDX, 
        :V_3_D_IDX, :V_3_Q_IDX, 
        :I_1_D_IDX, :I_1_Q_IDX, 
        :I_3_D_IDX, :I_3_Q_IDX, 
        :I_B1_D_IDX, :I_B1_Q_IDX, 
        :I_B2_D_IDX, :I_B2_Q_IDX, 
        :I_B3_D_IDX, :I_B3_Q_IDX
    ]
    
    # Machine constants
    machine_const = [
        :DELTA, :OMEGA, 
        :EQ_P, :ED_P, 
        :PSI_D_PP, :PSI_Q_PP, 
        :PSI_D, :PSI_Q
    ]
    
    # AVR constants
    avr_const = [
        :EFD_IDX, :VS_IDX, 
        :VLL_IDX, :VF_IDX
    ]
    
    # Governor constants
    governor_const = [
        :FV_IDX, :FF_IDX, :ET_IDX
    ]
    
    # Calculate number of states for each component
    num_network = 22
    num_machine = 8
    num_avr = 4
    num_governor = 3
    
    # Calculate starting indices
    network_start = 1
    machine_start = network_start + num_network
    avr_start = machine_start + num_machine
    governor_start = avr_start + num_avr
    
    # Create ranges
    network_range = network_start:(network_start + num_network - 1)
    machine_range = machine_start:(machine_start + num_machine - 1)
    avr_range = avr_start:(avr_start + num_avr - 1)
    governor_range = governor_start:(governor_start + num_governor - 1)
    
    # Create mappings
    network_map = Dict(key => network_start + idx - 1 for (idx, key) in enumerate(network_const))
    machine_map = Dict(key => machine_start + idx - 1 for (idx, key) in enumerate(machine_const))
    avr_map = Dict(key => avr_start + idx - 1 for (idx, key) in enumerate(avr_const))
    governor_map = Dict(key => governor_start + idx - 1 for (idx, key) in enumerate(governor_const))
    
    # Calculate total number of states
    num_states = network_start + num_network + num_machine + num_avr + num_governor - 1
    
    return StateIndices(
        network_map, machine_map, avr_map, governor_map,
        network_range, machine_range, avr_range, governor_range,
        num_states
    )
end

end # module