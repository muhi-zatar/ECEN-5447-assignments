using SynchronousMachineModels
using Plots

"""
Example demonstrating how to create a custom power system with modified parameters.
"""
function custom_system_example()
    # Create custom machine parameters
    machine_params = MachineParams(
        4.0,    # H - Higher inertia
        1.5,    # D - Lower damping
        1.9,    # Xd
        1.8,    # Xq
        0.32,   # Xdp
        0.55,   # Xqp
        0.25,   # Xdpp
        0.25,   # Xqpp
        0.15,   # Xl
        9.0,    # Td0p
        0.5,    # Tq0p
        0.03,   # Td0pp
        0.05,   # Tq0pp
        0.003   # Ra
    )
    
    # Create custom AVR parameters (higher gain)
    avr_params = AVRParams(
        300.0,  # Ka - Higher gain
        0.015,  # Ta - Faster response
        1.0,    # Ke
        0.2,    # Te
        0.05,   # Kf
        1.0,    # Tf
        6.0,    # Vr_max - Higher limit
        -6.0    # Vr_min - Lower limit
    )
    
    # Create custom turbine-governor parameters
    turbine_gov_params = TurbineGovParams(
        0.04,   # R - Lower droop (more sensitive)
        0.15,   # Tg - Faster governor
        0.25,   # T_ch - Faster steam chest
        6.0,    # T_rh - Faster reheat
        0.4,    # F_hp - Different power distribution
        0.6,    # F_lp
        1.2,    # P_max - Higher power limit
        0.0     # P_min
    )
    
    # Create custom network parameters (stronger connection)
    network_params = NetworkParams(
        0.0,    # R_e
        0.3,    # X_e - Stronger connection
        1.0     # V_∞
    )
    
    # Create components with custom parameters
    machine = SynchronousMachine(machine_params)
    avr = AVR(avr_params)
    turbine_gov = TurbineGovernor(turbine_gov_params)
    network = Network(network_params)
    
    # Create complete power system
    system = PowerSystem(machine, avr, turbine_gov, network)
    
    # Initialize system
    initialize!(system, P0=0.9)  # Higher initial power
    
    # Run simulations
    println("Running Custom System Voltage Step Test...")
    result_voltage = voltage_reference_step(system, 0.1, 1.0, 20.0)  # Larger step (10%)
    plot_voltage, _ = result_voltage
    savefig(plot_voltage, "custom_voltage_response.png")
    
    println("Running Custom System Load Change Test...")
    result_load = load_change_step(system, -0.15, 1.0, 25.0)  # Larger load change (15%)
    plot_load, _ = result_load
    savefig(plot_load, "custom_load_response.png")
    
    return system, plot_voltage, plot_load
end

# Run the custom system example
custom_system_example()