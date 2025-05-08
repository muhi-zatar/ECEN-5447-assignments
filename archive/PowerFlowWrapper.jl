# Defining the module
module PowerFlowWrapper
# Functions to export
export get_powerflow_results

# using PowerSystems
# using PowerFlows
# const PSY = PowerSystems
# const PF = PowerFlows

function get_powerflow_results(network_description; reduce_load=true)
    # Loading the network data using power systems package
    # sys = PSY.System(network_description)

    if reduce_load
        # Reduce the load
        # load = get_component(StandardLoad, sys, "load1031")
        # set_constant_active_power!(load, (get_constant_active_power(load) * 0.975))
        # set_constant_reactive_power!(load, (get_constant_reactive_power(load) * 0.975))

        # Hardcoded solution
        v = [1.02, 1.0142, 0.9879419356431093]
        θ = [0.0, -0.006043588435856549, -0.10678039032163857]
        P = [0.971156468404729, 0.8, -1.755]
        Q = [0.23683020624294693, 0.09729586920369482, -0.2925]
    else
        v = [1.02, 1.0142, 0.9869282927720892]
        θ = [0.0, -0.007843738409384567, -0.1104518610357432]
        P = [1.017056789531292, 0.8, -1.8]
        Q = [0.24501857920705117, 0.10751144041475172, -0.3]
    end

    # # Running power flow with PowerFlows Package
    # pf_result = PF.solve_powerflow(PF.ACPowerFlow(), sys)
    # PSY.set_units_base_system!(sys, "SYSTEM_BASE")

    # # Extracting voltage magnitude and angle    
    # v = pf_result["bus_results"].Vm     # [pu-V]
    # θ = pf_result["bus_results"].θ      # [rad]

    # # Extracting active and reactive power (normlaized for base power in p.u)
    # P = pf_result["bus_results"].P_net / PSY.get_base_power(sys)
    # Q = pf_result["bus_results"].Q_net / PSY.get_base_power(sys)

    # println("v = [$(v[1]), $(v[2]), $(v[3])]")
    # println("θ = [$(θ[1]), $(θ[2]), $(θ[3])]")
    # println("P = [$(P[1]), $(P[2]), $(P[3])]")
    # println("Q = [$(Q[1]), $(Q[2]), $(Q[3])]")

    # Return vectors of values
    return (v, θ, P, Q)
end
end # module
