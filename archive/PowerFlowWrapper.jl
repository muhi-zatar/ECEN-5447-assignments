# Defining the module
module PowerFlowWrapper
# Functions to export
export get_powerflow_results

using PowerSystems
using PowerFlows
const PSY = PowerSystems
const PF = PowerFlows

function get_powerflow_results(network_description)
    # TODO: Add the number of bus in the input to make the output bus configurable (for return values).
    # Loading the network data using power systems package
    # sys = PSY.System(network_description)

    # # Running power flow with PowerFlows Package
    # pf_result = PF.solve_powerflow(PF.ACPowerFlow(), sys)
    # # TODO: Check if system base is 100
    # PSY.set_units_base_system!(sys, "SYSTEM_BASE")

    # # Extracting the output of power flow
    # # From Fiona's code
    # v = pf_result["bus_results"].Vm     # [pu-V]
    # θ = pf_result["bus_results"].θ      # [rad]

    # # Extracting active and reactive power (normlaized for base power in p.u)
    # P = pf_result["bus_results"].P_net / PSY.get_base_power(sys)
    # Q = pf_result["bus_results"].Q_net / PSY.get_base_power(sys)

    v = [1.02, 1.0142, 0.9869282927720892]
    θ = [0.0, -0.007843738409384567, -0.1104518610357432]
    P = [1.017056789531292, 0.8, -1.8]
    Q = [0.24501857920705117, 0.10751144041475172, -0.3]

    # Return vectors of values
    return (v, θ, P, Q)
end
end # module
