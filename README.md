# ECEN-5447-assignments

PowerSystemDynamics/
├── Project.toml          # Julia package manager file
├── config/               # Configuration files
│   ├── default.toml      # Default configuration
│   └── scenarios/        # Different test scenarios
│       ├── line_trip.toml
│       ├── load_decrease.toml
│       └── load_increase.toml
├── data/                 # Power system data files
│   └── ThreeBusMultiLoad.raw
├── logs/                 # Log files directory
├── results/              # Simulation results and plots
├── scripts/              # Utility scripts
│   └── run_simulation.jl # Main execution script
└── src/                  # Source code
    ├── PowerSystemDynamics.jl      # Main module file
    ├── components/                 # Component models
    │   ├── AVR.jl                  # Automatic Voltage Regulator models
    │   ├── Governor.jl             # Governor models
    │   ├── Machine.jl              # Machine models
    │   └── Network.jl              # Network models
    ├── simulation/                 # Simulation related code
    │   ├── ComponentIndex.jl       # Component indexing
    │   ├── Config.jl               # Configuration handling
    │   ├── Initialization.jl       # System initialization
    │   ├── Perturbation.jl         # System perturbations
    │   └── Simulation.jl           # Simulation engine
    ├── utils/                      # Utility functions
    │   ├── Logging.jl              # Logging utilities
    │   ├── Plotting.jl             # Plotting utilities
    │   └── PowerFlow.jl            # Power flow wrapper
    └── types/                      # Type definitions
        └── SystemTypes.jl          # System type definitions