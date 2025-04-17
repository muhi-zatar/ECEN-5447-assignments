## References
- The model will draw from the Yazdani Voltage Source Converters text and the PowerSimulationsDynamics library

## Components
The inverter model will include the following components:
- An LCL filter to interface with the network, based on [these equations](https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/component_models/filters/#LCL-Filter-[LCLFilter])
- A reduced Order PLL based on [these equations](https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/component_models/freq_esti/)
- An outer loop controller for active and reactive power based on [these equations](https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/component_models/outer_control/#Active-Power-Droop-(P-droop)-and-Q-droop-[OuterControl])
- An innter loop controller for voltage control based on [these equations](https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/component_models/inner_control/#Integrated-Virtual-Impedance,-Voltage-and-Current-Controller-[VoltageModeControl])