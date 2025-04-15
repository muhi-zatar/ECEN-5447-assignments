# Main module for the inverter model
module InverterModel

# Exporting all component models
export FilterModel, PLLModel, DCLinkModel, OuterLoopModel, InnerLoopModel, PWMModel
export initialize_filter, update_filter_states!
export initialize_pll, update_pll_states!
export initialize_dclink, update_dclink_states!
export initialize_outerloop, update_outerloop_states!
export initialize_innerloop, update_innerloop_states!
export initialize_pwm, update_pwm_states!

# Include all component files
include("FilterModel.jl")
include("PLLModel.jl")
include("DCLinkModel.jl")
include("OuterLoopModel.jl")
include("InnerLoopModel.jl")
include("PWMModel.jl")

# Re-export components
using .FilterModel
using .PLLModel
using .DCLinkModel
using .OuterLoopModel
using .InnerLoopModel
using .PWMModel

end # module InverterModel