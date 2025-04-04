#!/usr/bin/env julia

# Run this script from the project root directory
# Usage: julia scripts/run_simulation.jl [config_file] [network_file]

using PowerSystemDynamics
using PowerSystemDynamics.ComponentIndex: create_state_indices
using PowerSystemDynamics.Plotting: plot_results
using PowerSystemDynamics.MachineComponents
using PowerSystemDynamics.AVRComponents
using PowerSystemDynamics.GovernorComponents
using PowerSystemDynamics.NetworkComponents
using PowerSystemDynamics.PowerFlowComponents
using PowerSystemDynamics.Logging: init_logger, log_info, close_logger
using PowerSystemDynamics.Config: load_config
using PowerSystemDynamics.Simulation: initialize_system, run_simulation
using ArgParse

function parse_commandline()
    s = ArgParseSettings(
        description = "Run power system dynamics simulation",
        version = "1.0",
        add_version = true
    )
    
    @add_arg_table! s begin
        "config"
            help = "Configuration file path"
            arg_type = String
            default = "config/default.toml"
        "network"
            help = "Network data file path"
            arg_type = String
            default = "data/ThreeBusMultiLoad.raw"
        "--no-plots", "-n"
            help = "Disable plot generation"
            action = :store_true
        "--output-dir", "-o"
            help = "Output directory for results and plots"
            arg_type = String
            default = "results"
        "--log-dir", "-l"
            help = "Directory for log files"
            arg_type = String
            default = "logs"
        "--disable-avr"
            help = "Disable the AVR"
            action = :store_true
        "--disable-governor"
            help = "Disable the governor"
            action = :store_true
        "--verbose", "-v"
            help = "Enable verbose output"
            action = :store_true
    end
    
    return parse_args(s)
end

function main()
    # Parse command line arguments
    args = parse_commandline()
    
    # Initialize logger
    logger = init_logger(args["log-dir"])    
    log_info(logger, "Starting simulation...")
    log_info(logger, "Using configuration file: $(args["config"])")
    log_info(logger, "Using network file: $(args["network"])")
    
    # Load configuration
    model, config = load_config(args["config"])
    
    # Apply command line overrides
    if args["disable-avr"]
        model.enabled_components["avr"] = false
        log_info(logger, "AVR disabled by command line argument")
    end
    
    if args["disable-governor"]
        model.enabled_components["governor"] = false
        log_info(logger, "Governor disabled by command line argument")
    end
    
    # Create component instances based on configuration
    model.network = create_network_from_config(config)
    model.machine = create_machine_from_config(config)
    model.avr = create_avr_from_config(config)
    model.governor = create_governor_from_config(config)
    
    # Initialize the system
    # log_info(logger, "Initializing system...")
    # initialize_system(model, args["network"], logger)
    
    # Run simulation
    log_info(logger, "Running simulation...")
    sol = run_simulation(model, logger)
    
    # Generate plots
    if !args["no-plots"]
        log_info(logger, "Generating plots...")
        indices = create_state_indices(model)
        plot_results(sol, indices, args["output-dir"])
    end
    
    # Close logger
    log_info(logger, "Simulation completed successfully")
    close_logger(logger)
    
    return sol
end

# Create components from config
function create_network_from_config(config)
    network_config = get(config, "network", Dict())
    
    return ThreeBusNetwork(
        R_12 = get(network_config, "R_12", 0.01),
        X_12 = get(network_config, "X_12", 0.12),
        B_1 = get(network_config, "B_1", 0.05),
        R_13 = get(network_config, "R_13", 0.01),
        X_13 = get(network_config, "X_13", 0.12),
        B_3 = get(network_config, "B_3", 0.05),
        R_23 = get(network_config, "R_23", 0.01),
        X_23 = get(network_config, "X_23", 0.12),
        B_2 = get(network_config, "B_2", 0.05),
        X_IB = get(network_config, "X_IB", 0.1)
    )
end

function create_machine_from_config(config)
    machine_config = get(config, "machine", Dict())
    
    return SauerPaiMachine(
        R = get(machine_config, "R", 0.002),
        X_d = get(machine_config, "X_d", 1.79),
        X_q = get(machine_config, "X_q", 1.71),
        Xd_p = get(machine_config, "Xd_p", 0.169),
        Xq_p = get(machine_config, "Xq_p", 0.228),
        Xd_pp = get(machine_config, "Xd_pp", 0.135),
        Xq_pp = get(machine_config, "Xq_pp", 0.2),
        Xl = get(machine_config, "Xl", 0.13),
        Td0_p = get(machine_config, "Td0_p", 4.3),
        Tq0_p = get(machine_config, "Tq0_p", 0.85),
        Td0_pp = get(machine_config, "Td0_pp", 0.032),
        Tq0_pp = get(machine_config, "Tq0_pp", 0.05),
        H = get(machine_config, "H", 3.148),
        D = get(machine_config, "D", 2.0),
        base_power = get(machine_config, "base_power", 100.0),
        system_base_power = get(machine_config, "system_base_power", 100.0),
        system_base_frequency = get(machine_config, "system_base_frequency", 60.0)
    )
end

function create_avr_from_config(config)
    avr_config = get(config, "avr", Dict())
    
    return EXST1(
        TR = get(avr_config, "TR", 0.01),
        TB = get(avr_config, "TB", 20.0),
        TC = get(avr_config, "TC", 10.0),
        KF = get(avr_config, "KF", 0.0),
        TF = get(avr_config, "TF", 0.1),
        KA = get(avr_config, "KA", 200.0),
        TA = get(avr_config, "TA", 0.1),
        KC = get(avr_config, "KC", 0.0),
        V_ref = get(avr_config, "V_ref", 1.0)
    )
end

function create_governor_from_config(config)
    gov_config = get(config, "governor", Dict())
    
    return GasTG(
        R = get(gov_config, "R", 0.05),
        T1 = get(gov_config, "T1", 0.2),
        T2 = get(gov_config, "T2", 0.2),
        T3 = get(gov_config, "T3", 2.0),
        D_turb = get(gov_config, "D_turb", 0.0),
        AT = get(gov_config, "AT", 1.0),
        KT = get(gov_config, "KT", 2.5),
        V_min = get(gov_config, "V_min", 0.01),
        V_max = get(gov_config, "V_max", 1.1),
        P_ref = get(gov_config, "P_ref", 0.8)
    )
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end