# PowerSystemDynamics

## Overview

PowerSystemDynamics is a Julia package for simulating power system dynamics. The package currently focuses on a three-bus power system with:

- An infinite bus (Bus 1)
- A generator bus with full dynamic models (Bus 2)
- A load bus (Bus 3)

The implementation includes detailed models for:
- Synchronous machine (Sauer-Pai model)
- Automatic Voltage Regulator (EXST1 type)
- Governor (Gas Turbine model)
- Network dynamics and algebraic constraints

## Project Structure

```
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
    │   ├── AVR.jl                  # AVR models (AVRComponents)
    │   ├── Governor.jl             # Governor models (GovernorComponents)
    │   ├── Machine.jl              # Machine models (MachineComponents)
    │   └── Network.jl              # Network models (NetworkComponents)
    ├── simulation/                 # Simulation related code
    │   ├── ComponentIndex.jl       # Component indexing
    │   ├── Config.jl               # Configuration handling
    │   ├── Initialization.jl       # System initialization
    │   ├── Perturbation.jl         # System perturbations
    │   └── Simulation.jl           # Simulation engine
    ├── utils/                      # Utility functions
    │   ├── Logging.jl              # Logging utilities
    │   ├── Plotting.jl             # Plotting utilities
    │   ├── PowerFlow.jl            # Power flow wrapper (PowerFlowComponents)
    │   └── Transformations.jl      # Coordinate transformation functions
    └── types/                      # Type definitions
        └── SystemTypes.jl          # System type definitions
```

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd PowerSystemDynamics
```

2. Activate and instantiate the Julia project:
```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Running Simulations

The package includes configuration files for common scenarios such as line trips and load changes.

To run a simulation:

```bash
julia --project=. scripts/run_simulation.jl config/scenarios/line_trip.toml
```

### Available Simulation Scenarios

1. **Line Trip**: Simulates the disconnection of a transmission line
   ```bash
   julia --project=. scripts/run_simulation.jl config/scenarios/line_trip.toml
   ```

2. **Load Decrease**: Simulates a decrease in load demand
   ```bash
   julia --project=. scripts/run_simulation.jl config/scenarios/load_decrease.toml
   ```

3. **Load Increase**: Simulates an increase in load demand
   ```bash
   julia --project=. scripts/run_simulation.jl config/scenarios/load_increase.toml
   ```

### Command Line Options

The `run_simulation.jl` script supports several command line options:

```
Usage: run_simulation.jl [options] config network

Arguments:
  config                    Configuration file path (default: "config/default.toml")
  network                   Network data file path (default: "data/ThreeBusMultiLoad.raw")

Options:
  -n, --no-plots            Disable plot generation
  -o, --output-dir DIR      Output directory for results and plots (default: "results")
  -l, --log-dir DIR         Directory for log files (default: "logs")
  --disable-avr             Disable the AVR
  --disable-governor        Disable the governor
  -v, --verbose             Enable verbose output
  -h, --help                Show this help
```

## Configuration

The simulation is controlled by TOML configuration files located in the `config/` directory. You can customize parameters for:

- Simulation settings (time span, step size, etc.)
- Component enablement (turn on/off AVR, governor, etc.)
- Network parameters (line impedances, shunt elements)
- Machine parameters (reactances, time constants)
- AVR parameters (gains, time constants)
- Governor parameters (droop, time constants)
- Perturbation settings (type, magnitude, timing)

Example configuration file:

```toml
# Simulation parameters
[simulation]
start_time = 0.0
end_time = 20.0
time_step = 0.001
perturb_time = 5.0
save_interval = 0.01

# Component enablement
[enabled_components]
machine = true
avr = true
governor = true

# Perturbation settings
[perturbation]
type = "LINE_TRIP"
line = "1-2"
```

## Development

### Adding New Components

To add new components to the system:

1. Create a new file in the appropriate directory (e.g., `src/components/`)
2. Define a new module with the component logic
3. Update the main module (`PowerSystemDynamics.jl`) to include and export the new component
4. Update the simulation logic to use the new component

### Testing

Basic tests can be run with:

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

## Output

Simulation results are saved to:

- **Logs**: `logs/` directory contains detailed simulation logs
- **Plots**: `results/` directory contains generated plots

## Acknowledgments

This package builds on work by several open-source projects and academic resources, including:
- DifferentialEquations.jl for ODE/DAE solving
- PowerSystems.jl and PowerFlows.jl for power flow solutions (with adapters)
- Milano's "Power System Modelling and Scripting" for theoretical foundations
