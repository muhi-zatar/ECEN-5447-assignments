module PowerSystemDynamics

# Export all public interfaces
export run_simulation, load_config, plot_results, initialize_system
export SauerPaiMachine, EXST1, GasTG, ThreeBusNetwork
export apply_perturbation!, Logger
export ri_dq, dq_ri
export BUS_INFINITE_BUS, BUS_MACHINE_MODEL, BUS_LOAD

# Core types
include("types/SystemTypes.jl")
include("utils/Transformations.jl")

# Utility transformation functions
include("components/Machine.jl")    # MachineComponents
include("components/AVR.jl")        # AVRComponents
include("components/Governor.jl")   # GovernorComponents
include("components/Network.jl")    
include("simulation/ComponentIndex.jl")
include("utils/Logging.jl")
include("utils/Plotting.jl")
include("utils/PowerFlow.jl")
# Component models - with explicit module names to avoid conflicts
# NetworkComponents

# Simulation utilities
include("simulation/Config.jl")
include("simulation/Initialization.jl")
include("simulation/Perturbation.jl")
include("simulation/Simulation.jl")

# Utility functions

# Import and re-export from submodules
import .Types
import .Transformations
import .ComponentIndex
import .PowerFlowComponents
import .MachineComponents: SauerPaiMachine, initialize_machine_states, update_machine_states!
import .MachineComponents: DELTA, OMEGA, EQ_P, ED_P, PSI_D_PP, PSI_Q_PP, PSI_D, PSI_Q
import .AVRComponents: EXST1, initialize_avr_states, update_avr_states!
import .AVRComponents: EFD_IDX, VS_IDX, VLL_IDX, VF_IDX
import .GovernorComponents: GasTG, initialize_gov_states, update_gov_states!
import .GovernorComponents: FV_IDX, FF_IDX, ET_IDX
import .NetworkComponents: ThreeBusNetwork, initialize_network_states, update_network_states!
import .NetworkComponents: BUS_INFINITE_BUS, BUS_MACHINE_MODEL, BUS_LOAD
import .NetworkComponents: I_12_D_IDX, I_12_Q_IDX, I_13_D_IDX, I_13_Q_IDX, I_23_D_IDX, I_23_Q_IDX
import .NetworkComponents: V_1_D_IDX, V_1_Q_IDX, V_2_D_IDX, V_2_Q_IDX, V_3_D_IDX, V_3_Q_IDX
import .NetworkComponents: I_1_D_IDX, I_1_Q_IDX, I_3_D_IDX, I_3_Q_IDX
import .NetworkComponents: I_B1_D_IDX, I_B1_Q_IDX, I_B2_D_IDX, I_B2_Q_IDX, I_B3_D_IDX, I_B3_Q_IDX
import .NetworkComponents: sanity_check
import .PowerFlowComponents: get_powerflow_results
import .Transformations: ri_dq, dq_ri
import .Config
import .Initialization
import .Perturbation
import .Simulation
import .Logging
import .Plotting

# Make component modules directly accessible
const Machine = MachineComponents
const AVR = AVRComponents
const Governor = GovernorComponents
const Network = NetworkComponents
const PowerFlow = PowerFlowComponents

end # module