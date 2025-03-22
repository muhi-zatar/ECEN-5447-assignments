using PowerSystemCaseBuilder
using PowerSimulationsDynamics
using PowerNetworkMatrices
using Plots
using Sundials
using PowerSystems

sys = build_system(PSIDTestSystems, "psid_test_threebus_multimachine_dynlines")

tspan = (0.0, 30.0)

sim = Simulation(
           ResidualModel, #Type of model used
           sys, #system
           pwd(), #folder to output results
           tspan, #time span
           # Ybus_change, #Type of perturbation
       )

execute!(
    sim, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = 0.02, #Maximum step size
    )

results = read_results(sim)