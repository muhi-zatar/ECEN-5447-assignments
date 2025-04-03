module Types

export PowerSystemModel, SimulationParams, BusType

# Enum for bus types
@enum BusType begin
    INFINITE_BUS = 1
    MACHINE_BUS = 2
    LOAD_BUS = 3
end

"""
    SimulationParams

Contains parameters for the simulation run.
"""
mutable struct SimulationParams
    tspan::Tuple{Float64, Float64}  # Time span for simulation (start, end)
    dt::Float64                     # Time step
    perturb_time::Float64           # Time to apply perturbation
    save_interval::Float64          # Interval to save results
    
    # Constructor with default values
    function SimulationParams(;
        tspan=(0.0, 30.0),
        dt=0.0001,
        perturb_time=5.0,
        save_interval=0.01
    )
        return new(tspan, dt, perturb_time, save_interval)
    end
end

"""
    PowerSystemModel

Container for all the components of the power system.
"""
mutable struct PowerSystemModel
    network::Any                    # Network model
    machine::Any                    # Machine model
    avr::Any                        # AVR model
    governor::Any                   # Governor model
    enabled_components::Dict{String, Bool}  # Enabled components
    simulation_params::SimulationParams     # Simulation parameters
    
    # Constructor with default values
    function PowerSystemModel(;
        network=nothing,
        machine=nothing,
        avr=nothing,
        governor=nothing,
        enabled_components=Dict(
            "machine" => true,
            "avr" => true,
            "governor" => true
        ),
        simulation_params=SimulationParams()
    )
        return new(network, machine, avr, governor, enabled_components, simulation_params)
    end
end

end # module