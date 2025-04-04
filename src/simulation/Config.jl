module Config

export load_config, save_config

using TOML
using ..Types

"""
    load_config(filename::String)

Load configuration from a TOML file and create a PowerSystemModel.
"""
function load_config(filename::String)
    # Check if file exists
    if !isfile(filename)
        error("Configuration file not found: $filename")
    end
    
    # Parse TOML file
    config = TOML.parsefile(filename)
    
    # Create simulation parameters
    sim_params = get_simulation_params(config)
    
    # Component enablement
    enabled = get(config, "enabled_components", Dict{String, Bool}())
    enabled_components = Dict(
        "machine" => get(enabled, "machine", true),
        "avr" => get(enabled, "avr", true),
        "governor" => get(enabled, "governor", true)
    )
    
    # Create model shell
    model = PowerSystemModel(
        enabled_components=enabled_components,
        simulation_params=sim_params
    )
    
    return model, config
end

"""
    get_simulation_params(config::Dict)

Extract simulation parameters from config dictionary.
"""
function get_simulation_params(config::Dict)
    # Get simulation section or use empty dict if not present
    sim_config = get(config, "simulation", Dict())
    
    # Extract parameters with defaults
    tspan = (
        get(sim_config, "start_time", 0.0),
        get(sim_config, "end_time", 20.0)
    )
    dt = get(sim_config, "time_step", 0.001)
    perturb_time = get(sim_config, "perturb_time", 5.0)
    save_interval = get(sim_config, "save_interval", 0.01)
    
    return SimulationParams(
        tspan=tspan,
        dt=dt,
        perturb_time=perturb_time,
        save_interval=save_interval
    )
end

"""
    save_config(filename::String, model::PowerSystemModel, additional_data::Dict=Dict())

Save the current configuration to a TOML file.
"""
function save_config(filename::String, model::PowerSystemModel, additional_data::Dict=Dict())
    # Create config dictionary
    config = Dict{String, Any}()
    
    # Add enabled components
    config["enabled_components"] = model.enabled_components
    
    # Add simulation parameters
    config["simulation"] = Dict(
        "start_time" => model.simulation_params.tspan[1],
        "end_time" => model.simulation_params.tspan[2],
        "time_step" => model.simulation_params.dt,
        "perturb_time" => model.simulation_params.perturb_time,
        "save_interval" => model.simulation_params.save_interval
    )
    
    # Add additional data
    for (key, value) in additional_data
        config[key] = value
    end
    
    # Create directory if it doesn't exist
    dir = dirname(filename)
    if !isdir(dir)
        mkpath(dir)
    end
    
    # Write to file
    open(filename, "w") do io
        TOML.print(io, config)
    end
    
    return config
end

end # module