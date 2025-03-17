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
    sys = PSY.System(network_description)

    # Running power flow with PowerFlows Package
    pf_result = PF.solve_powerflow(PF.ACPowerFlow(), sys)
    # TODO: Check if system base is 100
    PSY.set_units_base_system!(sys, "SYSTEM_BASE")

    # Extracting the output of power flow
    # From Fiona's code
    v = pf_result["bus_results"].Vm     # [pu-V]
    θ = pf_result["bus_results"].θ      # [rad]

    # Extracting active and reactive power (normlaized for base power in p.u)
    P = pf_result["bus_results"].P_net / PSY.get_base_power(sys)
    Q = pf_result["bus_results"].Q_net / PSY.get_base_power(sys)

    # Return values for Bus 2 (where we are placing our model :D)
    return (v[2], θ[2], P[2], Q[2])
end
end # module