# Defining the module
module PerturbationModel

# Exporting the necessary functions
export Perturbation, apply_perturbation!

# Define perturbation types
const LOAD_CHANGE = 1
const FAULT = 2
const REFERENCE_CHANGE = 3
const LINE_TRIP = 4

# Perturbation structure
mutable struct Perturbation
    type::String
    start_time::Float64
    end_time::Float64
    bus_id::Int
    branch_from::Int
    branch_to::Int
    magnitude::Float64
    
    # Constructor for load changes and reference changes
    function Perturbation(type, start_time, end_time, bus_id, magnitude)
        return new(type, start_time, end_time, bus_id, 0, 0, magnitude)
    end
end

# Function to apply a perturbation to the system
function apply_perturbation!(network, perturbation, current_time)
    # Check if perturbation is active
    is_active = current_time >= perturbation.start_time && current_time < perturbation.end_time
    
    if perturbation_load.type == "LOAD_CHANGE" && is_active
        # Apply load change to specified bus
        bus = network.buses[perturbation.bus_id]
        # Store original load if this is the first time applying
        if current_time == perturbation.start_time
            perturbation.original_load_P = bus.P_load
            perturbation.original_load_Q = bus.Q_load
        end
        
        # Change load by specified magnitude
        bus.P_load = perturbation.original_load_P * (1.0 + perturbation.magnitude)
        bus.Q_load = perturbation.original_load_Q * (1.0 + perturbation.magnitude)
        
    elseif perturbation.type == "LOAD_CHANGE" && current_time == perturbation.end_time
        # Restore original load
        bus = network.buses[perturbation.bus_id]
        bus.P_load = perturbation.original_load_P
        bus.Q_load = perturbation.original_load_Q
        
    elseif perturbation.type == REFERENCE_CHANGE && is_active
        # Apply reference voltage/power change to a generator
        if perturbation.bus_id > 0
            bus = network.buses[perturbation.bus_id]
            
            # Store original values if this is the first time applying
            if current_time == perturbation.start_time
                if bus.type == PV
                    perturbation.original_value = bus.V_mag
                else
                    perturbation.original_value = bus.P_gen
                end
            end
            
            # Apply change
            if bus.type == PV
                bus.V_mag = perturbation.original_value * (1.0 + perturbation.magnitude)
            else
                bus.P_gen = perturbation.original_value * (1.0 + perturbation.magnitude)
            end
        end
        
    elseif perturbation.type == REFERENCE_CHANGE && current_time == perturbation.end_time
        # Restore original reference
        bus = network.buses[perturbation.bus_id]
        if bus.type == PV
            bus.V_mag = perturbation.original_value
        else
            bus.P_gen = perturbation.original_value
        end
        
    elseif perturbation.type == FAULT && is_active
        # Apply fault (simplified as an impedance to ground)
        # This would require modifying the Y-bus matrix temporarily
        # TODO: modify the Y-bus if this fault is needed.
        println("Fault applied at time $current_time on branch $(perturbation.branch_from)-$(perturbation.branch_to)")
        
    elseif perturbation.type == LINE_TRIP && is_active
        # Trip a line by setting its admittance to zero
        # Find the branch
        for (i, branch) in enumerate(network.branches)
            if branch.from_bus == perturbation.branch_from && branch.to_bus == perturbation.branch_to
                # Store original values if this is the first time applying
                if current_time == perturbation.start_time
                    perturbation.original_r = branch.r
                    perturbation.original_x = branch.x
                end
                
                # Set impedance to very high values (effectively removing the branch)
                # TODO: Check if this is the way to do it.
                branch.r = 1e6
                branch.x = 1e6
                println("Line $(perturbation.branch_from)-$(perturbation.branch_to) tripped at time $current_time")
                break
            end
        end
        
    elseif perturbation.type == LINE_TRIP && current_time == perturbation.end_time
        # Restore line
        for (i, branch) in enumerate(network.branches)
            if branch.from_bus == perturbation.branch_from && branch.to_bus == perturbation.branch_to
                branch.r = perturbation.original_r
                branch.x = perturbation.original_x
                println("Line $(perturbation.branch_from)-$(perturbation.branch_to) restored at time $current_time")
                break
            end
        end
    end
    
    return is_active
end

end # module