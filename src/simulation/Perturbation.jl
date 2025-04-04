module Perturbation

export apply_perturbation!, PerturbationType

using ..Types

"""
    PerturbationType

Enum for different types of system perturbations.
"""
@enum PerturbationType begin
    LINE_TRIP       # Disconnect a transmission line
    LOAD_INCREASE   # Increase the load at a bus
    LOAD_DECREASE   # Decrease the load at a bus
end

"""
    apply_perturbation!(model::PowerSystemModel; 
                        perturbation_type::PerturbationType=LINE_TRIP,
                        parameter_changes::Dict{Symbol, Any}=Dict())

Apply a perturbation to the power system model.
"""
function apply_perturbation!(model::PowerSystemModel; 
                            perturbation_type::PerturbationType=LINE_TRIP,
                            parameter_changes::Dict{Symbol, Any}=Dict())
    # Apply perturbation based on type
    if perturbation_type == LINE_TRIP
        apply_line_trip!(model, parameter_changes)
    elseif perturbation_type == LOAD_INCREASE
        apply_load_change!(model, parameter_changes, true)
    elseif perturbation_type == LOAD_DECREASE
        apply_load_change!(model, parameter_changes, false)
    else
        @warn "Unknown perturbation type: $perturbation_type"
    end
end

"""
    apply_line_trip!(model::PowerSystemModel, params::Dict{Symbol, Any})

Apply a line trip perturbation to the model.
"""
function apply_line_trip!(model::PowerSystemModel, params::Dict{Symbol, Any})
    # Default to line 1-2 if not specified
    line = get(params, :line, "1-2")
    
    # Apply perturbation based on which line is tripped
    if line == "1-2"
        # Make line 1-2 very high impedance (effectively open circuit)
        model.network.R_12 = 1e6
        model.network.X_12 = 1e6
        # Update shunt capacitance to maintain numerical stability
        model.network.B_1 *= 0.5
        model.network.B_2 *= 0.5
    elseif line == "1-3"
        model.network.R_13 = 1e6
        model.network.X_13 = 1e6
        # Update shunt capacitance
        model.network.B_1 *= 0.5
        model.network.B_3 *= 0.5
    elseif line == "2-3"
        model.network.R_23 = 1e6
        model.network.X_23 = 1e6
        # Update shunt capacitance
        model.network.B_2 *= 0.5
        model.network.B_3 *= 0.5
    else
        @warn "Unknown line: $line. No perturbation applied."
    end
end

"""
    apply_load_change!(model::PowerSystemModel, params::Dict{Symbol, Any}, increase::Bool)

Apply a load change perturbation to the model.
"""
function apply_load_change!(model::PowerSystemModel, params::Dict{Symbol, Any}, increase::Bool)
    # Default to load at bus 3 if not specified
    bus = get(params, :bus, "3")
    # Default to 15% change if not specified
    percentage = get(params, :percentage, 0.15)
    
    if increase
        # For load increase, decrease the load impedance
        factor = 1.0 / (1.0 + percentage)
    else
        # For load decrease, increase the load impedance
        factor = 1.0 + percentage
    end
    
    # Apply perturbation based on which bus has load change
    if bus == "3"
        # Change load impedance at bus 3
        model.network.Z_L *= factor
    else
        @warn "Unknown load bus: $bus. No perturbation applied."
    end
end

end # module