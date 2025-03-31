cd(@__DIR__)
using Pkg
Pkg.activate("../.")
#Pkg.resolve() # make sure Manifest matches Project
#Pkg.instantiate() # install missing dependencies

# Importing necessary Modules
using PowerSystemCaseBuilder
using PowerSimulationsDynamics
using PowerSystems
using Sundials
using Plots
using Plots.PlotMeasures
const PSY = PowerSystems

# Load the system
omib_sys = build_system(PSIDSystems, "OMIB System")

# Build the simulation and initialize the problem
time_span = (0.0, 30.0)

perturbation_trip = BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1")

sim = Simulation(
    ResidualModel, # Type of formulation
    omib_sys, # System
    mktempdir(), # Output directory
    time_span,
    perturbation_trip)

show_states_initial_value(sim)                      # Show the initial states

get_component(Source, omib_sys, "InfBus")           # Show the infinite bus metadata

# Run the simulation
execute!(
    sim, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax=0.02, #Arguments: Maximum timestep allowed
);

# Explore the solution
results = read_results(sim)

voltage_IB = get_voltage_magnitude_series(results, 101)
voltage_SM = get_voltage_magnitude_series(results, 102)

plot(voltage_IB, label="Inf. Bus", title="Voltage Response to Line Trip", xlabel="time", ylabel="Voltage [pu]", linewidth=2)
plot!(voltage_SM, label="SM")