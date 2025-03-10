module PowerFlowWrapper

export get_powerflow_results, solve_powerflow_with_perturbation

using PowerSystems
using PowerFlows
const PSY = PowerSystems
const PF = PowerFlows

"""
    get_powerflow_results(network_description)

Loads a network from a raw file and runs power flow analysis.
Returns the voltage magnitude, angle, active power, and reactive power for Bus 2.

# Arguments
- `network_description`: Path to the raw file containing network description

# Returns
- Tuple of (voltage magnitude, voltage angle, active power, reactive power) at Bus 2
"""
function get_powerflow_results(network_description)
    # Load raw file into Sienna
    sys = PSY.System(network_description)

    # Run power flow
    pf_result = PF.solve_powerflow(PF.ACPowerFlow(), sys)
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")

    # Get power flow output
    v = pf_result["bus_results"].Vm     # [pu-V]
    θ = pf_result["bus_results"].θ      # [rad]

    # Sienna exports power injections in MW/MVar, so adjust by system base
    P = pf_result["bus_results"].P_net / PSY.get_base_power(sys) # [pu(MW)]
    Q = pf_result["bus_results"].Q_net / PSY.get_base_power(sys) # [pu(MVar)]

    # Return values for Bus 2, where we're placing our model
    return (v[2], θ[2], P[2], Q[2])
end

"""
    solve_powerflow_with_perturbation(network_description, perturbation_type, parameters)

Runs power flow with a specific perturbation applied to the system.

# Arguments
- `network_description`: Path to the raw file containing network description
- `perturbation_type`: Type of perturbation (e.g., "load_change", "line_trip")
- `parameters`: Dictionary of perturbation parameters

# Returns
- Tuple of (voltage magnitude, voltage angle, active power, reactive power) at Bus 2 after perturbation
"""
function solve_powerflow_with_perturbation(network_description, perturbation_type, parameters)
    # Load raw file into Sienna
    sys = PSY.System(network_description)
    
    # Apply perturbation based on type
    if perturbation_type == "load_change" && haskey(parameters, "bus_number") && haskey(parameters, "scaling_factor")
        bus_number = parameters["bus_number"]
        scaling_factor = parameters["scaling_factor"]
        
        # Find the load at the specified bus
        for load in PSY.get_components(PSY.PowerLoad, sys)
            if PSY.get_number(PSY.get_bus(load)) == bus_number
                # Scale the load
                original_p = PSY.get_active_power(load)
                original_q = PSY.get_reactive_power(load)
                PSY.set_active_power!(load, original_p * scaling_factor)
                PSY.set_reactive_power!(load, original_q * scaling_factor)
                break
            end
        end
    elseif perturbation_type == "line_trip" && haskey(parameters, "from_bus") && haskey(parameters, "to_bus")
        from_bus = parameters["from_bus"]
        to_bus = parameters["to_bus"]
        
        # Find the branch and set its status to false (tripped)
        for branch in PSY.get_components(PSY.Line, sys)
            from = PSY.get_number(PSY.get_from_bus(branch))
            to = PSY.get_number(PSY.get_to_bus(branch))
            
            if (from == from_bus && to == to_bus) || (from == to_bus && to == from_bus)
                PSY.set_available!(branch, false)
                break
            end
        end
    elseif perturbation_type == "voltage_change" && haskey(parameters, "bus_number") && haskey(parameters, "new_voltage")
        bus_number = parameters["bus_number"]
        new_voltage = parameters["new_voltage"]
        
        # Find the generator at the specified bus and change its voltage setpoint
        for gen in PSY.get_components(PSY.Generator, sys)
            if PSY.get_number(PSY.get_bus(gen)) == bus_number
                if hasmethod(PSY.set_voltage_setpoint!, Tuple{typeof(gen), Float64})
                    PSY.set_voltage_setpoint!(gen, new_voltage)
                    break
                end
            end
        end
    end
    
    # Run power flow with the perturbation applied
    pf_result = PF.solve_powerflow(PF.ACPowerFlow(), sys)
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")

    # Get power flow output
    v = pf_result["bus_results"].Vm     # [pu-V]
    θ = pf_result["bus_results"].θ      # [rad]

    # Sienna exports power injections in MW/MVar, so adjust by system base
    P = pf_result["bus_results"].P_net / PSY.get_base_power(sys) # [pu(MW)]
    Q = pf_result["bus_results"].Q_net / PSY.get_base_power(sys) # [pu(MVar)]

    # Return values for Bus 2
    return (v[2], θ[2], P[2], Q[2])
end

end # module