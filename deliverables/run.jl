# Ensure script runs in the correct environment
cd(@__DIR__)
using Pkg
Pkg.activate("../.")
#Pkg.resolve()
#Pkg.instantiate()

# Importing the required librarires
using DifferentialEquations
using Plots
using LinearAlgebra

# Import our custom modules from other files
include("NetworkModel.jl")
include("SauerPaiMachineModel.jl")
include("EXST1Model.jl")
include("GasTGModel.jl")
include("PowerFlowWrapper.jl")

# Use our modules
using .NetworkModel
using .SauerPaiMachineModel
using .EXST1Model
using .GasTGModel
using .PowerFlowWrapper


# Definining some constants
const BUS_INFINITE_BUS = 1
const BUS_MACHINE_MODEL = 2
const BUS_LOAD = 3
const NUM_STATES = 22
const I_12_D_IDX = 1
const I_12_Q_IDX = 2
const I_13_D_IDX = 3
const I_13_Q_IDX = 4
const I_23_D_IDX = 5
const I_23_Q_IDX = 6
const V_1_D_IDX = 7
const V_1_Q_IDX = 8
const V_2_D_IDX = 9
const V_2_Q_IDX = 10
const V_3_D_IDX = 11
const V_3_Q_IDX = 12
const I_1_D_IDX = 13
const I_1_Q_IDX = 14
const I_3_D_IDX = 15
const I_3_Q_IDX = 16
const I_B1_D_IDX = 17
const I_B1_Q_IDX = 18
const I_B2_D_IDX = 19
const I_B2_Q_IDX = 20
const I_B3_D_IDX = 21
const I_B3_Q_IDX = 22
const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6
const PSI_D = 7
const PSI_Q = 8
# AVR state indices
const EFD_IDX = 1   # Field voltage (output of amplifier)
const VS_IDX = 2    # Sensed terminal voltage
const VLL_IDX = 3   # Lead-lag output
const VF_IDX = 4    # Feedback signal
# Governor state indices
const FV_IDX = 1                # Fuel value
const FF_IDX = 2                # Fuel flow
const ET_IDX = 3                # Exhaust temp

# Matrix transformation functions
# These are used to convert between rectangular and polar coordinates
# Same as power system package to avoid confusion in this part.
function ri_dq(δ::T) where {T<:Number}
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return Float64[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

function dq_ri(δ::T) where {T<:Number}
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return Float64[
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end


# Define parameters for ODE solver
struct ODEParams
    network::ThreeBusNetwork
    machine::SauerPaiMachine
    avr::EXST1
    governor::GasTG
    network_idx::UnitRange{Int64}
    machine_idx::UnitRange{Int64}
    avr_idx::UnitRange{Int64}
    gov_idx::UnitRange{Int64}
end

# Main function to run the simulation
function run_simulation(network_file)
    # Step 1: Get initial conditions from power flow (all buses)
    V_sol, θ_sol, P_sol, Q_sol = get_powerflow_results(network_file)

    # Select bus for our model :D
    V_mag = V_sol[BUS_MACHINE_MODEL]
    V_angle = θ_sol[BUS_MACHINE_MODEL]
    P = P_sol[BUS_MACHINE_MODEL]
    Q = Q_sol[BUS_MACHINE_MODEL]

    # Convert voltage to complex form
    V_terminal_init = V_mag * exp(im * V_angle)
    I_terminal_init = conj(complex(P, Q) / V_terminal_init)
    println("Initial power flow results:")
    println("V_terminal = $V_mag p.u $V_angle rad")
    println("P = $P pu, Q = $Q pu")
    println("I_terminal = $(abs(I_terminal_init)) p.u., $(angle(I_terminal_init)) rad")

    # Step 2: Initialize models
    # Create model instances
    network = ThreeBusNetwork()

    machine = SauerPaiMachine(
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

    avr = EXST1(
        V_ref=V_mag  # Set reference voltage to initial voltage
    )

    governor = GasTG(
        P_ref=P  # Set reference power to initial power
    )

    # Initialize states
    network_states, _, _ = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)
    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal_init, V_angle, P, Q)
    avr_states = initialize_avr(avr, V_mag, Vf_init)  # Assume zero field current initially, is that right?
    ω_init = 1.0
    governor_states = initialize_gov(governor, τ_m_init, ω_init)

    # Combine all states into a single vector for ODE solver
    states = vcat(network_states, machine_states, avr_states, governor_states)

    # Print initial states for debugging
    println("\nInitial States:")
    println("Machine states: $machine_states")
    println("AVR states: $avr_states")
    println("Governor states: $governor_states")

    # Define state indices for easier access
    network_idx = 1:22
    machine_idx = 23:30
    avr_idx = 31:34
    gov_idx = 35:37

    # Create ODE parameters structure
    p = ODEParams(
        network,
        machine,
        avr,
        governor,
        network_idx,
        machine_idx,
        avr_idx,
        gov_idx,
    )

    # Build mass matrix
    M_system = zeros(Float64, length(states))

    if isa(network.M, Vector)
        M_system[network_idx] .= network.M
    else
        for i in 1:length(network_idx)
            M_system[network_idx[i]] = network.M[i, i]
        end
    end

    M_system[machine_idx] .= 1.0
    M_system[avr_idx] .= 1.0
    M_system[gov_idx] .= 1.0

    mass_matrix = Diagonal(M_system)
    println("Mass matrix = $(M_system)")

    # ODE function
    function ode_system!(du, u, p, t)
        # Extract states for each component
        network_states = u[p.network_idx]
        machine_states = u[p.machine_idx]
        avr_states = u[p.avr_idx]
        gov_states = u[p.gov_idx]

        # Extract parameters
        network = p.network
        machine = p.machine
        avr = p.avr
        governor = p.governor

        # Make copies of states for Float64 compatibility if needed
        network_states_f64 = convert.(Float64, network_states)
        machine_states_f64 = convert.(Float64, machine_states)
        avr_states_f64 = convert.(Float64, avr_states)
        gov_states_f64 = convert.(Float64, gov_states)

        # Arrays for derivatives
        du_network = zeros(Float64, length(p.network_idx))
        du_machine = zeros(Float64, length(p.machine_idx))
        du_avr = zeros(Float64, length(p.avr_idx))
        du_gov = zeros(Float64, length(p.gov_idx))

        # Grab terminal voltage from network
        v_2_d = network_states_f64[V_2_D_IDX]
        v_2_q = network_states_f64[V_2_Q_IDX]
        V_RI = dq_ri(0.0) * [v_2_d; v_2_q]
        V_terminal = complex(V_RI[1], V_RI[2])
        V_terminal_mag = abs(V_terminal)

        # Grab omega from machine
        ω = machine_states_f64[OMEGA]

        # Grab field voltage from avr
        efd = avr_states_f64[EFD_IDX]

        # Update the states of each component
        update_avr_states!(avr_states_f64, du_avr, V_terminal_mag, avr)

        τ_m = update_gov_states!(gov_states_f64, du_gov, ω, governor)

        _, S_terminal_machine, _, _, _ = update_machine_states!(machine_states_f64, du_machine, V_terminal, efd, τ_m, machine)

        _, _, _ = update_network_states!(network_states_f64, du_network, S_terminal_machine, network)

        # Copy derivatives to output
        du[p.network_idx] .= du_network
        du[p.machine_idx] .= du_machine
        du[p.avr_idx] .= du_avr
        du[p.gov_idx] .= du_gov

        # For debugging printing some results,
        if abs(t - round(t)) < 0.000001
            println("t=$t")
            println("Machine States: δ=$(machine_states_f64[DELTA]), ω=$(machine_states_f64[OMEGA]), EQ_P=$(machine_states_f64[EQ_P]), ED_P=$(machine_states_f64[ED_P]), PSI_D_PP=$(machine_states_f64[PSI_D_PP]), PSI_Q_PP=$(machine_states_f64[PSI_Q_PP])")
            println("AVR States: VF=$(avr_states_f64[EFD_IDX]), VT=$(avr_states_f64[VS_IDX]), VLL=$(avr_states_f64[VLL_IDX]), VFB=$(avr_states_f64[VF_IDX])")
            println("Gov States: FV=$(gov_states_f64[FV_IDX]), FF=$(gov_states_f64[FF_IDX]), ET=$(gov_states_f64[ET_IDX])")
            println("Network States: I_12_D=$(network_states_f64[I_12_D_IDX]), I_12_Q=$(network_states_f64[I_12_Q_IDX]), I_13_D=$(network_states_f64[I_13_D_IDX]), I_13_Q=$(network_states_f64[I_13_Q_IDX]), V_2_D=$(network_states_f64[V_2_D_IDX]), V_2_Q=$(network_states_f64[V_2_Q_IDX]), I_B2_D=$(network_states_f64[I_B2_D_IDX]), I_B2_Q=$(network_states_f64[I_B2_Q_IDX])")
        end
    end

    # Build the function
    explicitDAE_M = ODEFunction(ode_system!, mass_matrix=mass_matrix)

    # Step 4: Set up and solve the ODE system
    tspan = (0.0, 20.0)
    prob = ODEProblem(explicitDAE_M, states, tspan, p)

    # Define the set of times to apply a perturbation
    perturb_times = [5.0]

    # Define the condition for which to apply a perturbation
    function condition(u, t, integrator)
        t in perturb_times
    end

    # Define the perturbation
    function affect!(integrator)
        #### Uncomment the desired perturbation ####
        # Load Jump
        #integrator.p.network.Z_L *= 1.15

        # Load Decrease
        #integrator.p.network.Z_L *= 0.85

        # Line Trip
        integrator.p.network.R_12 = 1e6
        integrator.p.network.X_12 = 1e6
        integrator.p.network.B_1 *= 0.5
        integrator.p.network.B_2 *= 0.5
    end

    # Create a Callback function that represents the perturbation
    cb = DiscreteCallback(condition, affect!)

    # Run simulation
    sol = solve(prob, Rodas5(autodiff=false), dt=0.001, adaptive=false, saveat=0.01, callback=cb, tstops=perturb_times)

    # Process results
    t = sol.t
    delta_values = [sol[machine_idx[1], i] for i in 1:length(t)]
    omega_values = [sol[machine_idx[2], i] for i in 1:length(t)]
    Vf_values = [sol[avr_idx[1], i] for i in 1:length(t)]

    # Create plots - Just as an example, 90% of them are not needed
    # Plot shaft states
    p1 = plot(t, delta_values,
        label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, omega_values,
        label="Rotor speed (ω)", linewidth=2)
    savefig(p1, "../results/shaft_states.png")

    # Plot fluxes
    p2 = plot(t, [sol[machine_idx[3], i] for i in 1:length(t)],
        label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[machine_idx[4], i] for i in 1:length(t)],
        label="ed'", linewidth=2)
    plot!(p2, t, [sol[machine_idx[5], i] for i in 1:length(t)],
        label="ψd''", linewidth=2)
    plot!(p2, t, [sol[machine_idx[6], i] for i in 1:length(t)],
        label="ψq''", linewidth=2)
    savefig(p2, "../results/fluxes.png")

    # Plot avr states
    p3 = plot(t, Vf_values,
        label="Field Voltage", title="AVR States", linewidth=2)
    plot!(p3, t, [sol[avr_idx[2], i] for i in 1:length(t)],
        label="Terminal Voltage", linewidth=2)
    plot!(p3, t, [sol[avr_idx[3], i] for i in 1:length(t)],
        label="Lead-Lag State", linewidth=2)
    plot!(p3, t, [sol[avr_idx[4], i] for i in 1:length(t)],
        label="Feedback State", linewidth=2)
    savefig(p3, "../results/avr_states.png")

    # Plot governor states
    p4 = plot(t, [sol[gov_idx[1], i] for i in 1:length(t)],
        label="Fuel Value", title="Governor States", linewidth=2)
    plot!(p4, t, [sol[gov_idx[2], i] for i in 1:length(t)],
        label="Fuel Flow", linewidth=2)
    plot!(p4, t, [sol[gov_idx[3], i] for i in 1:length(t)],
        label="Exhaust Temp", linewidth=2)
    savefig(p4, "../results/gov_states.png")

    # Plot line currents
    p5 = plot(t, [sol[network_idx[I_12_D_IDX], i] for i in 1:length(t)], title="Line Currents", label="I_12_D", linewith=2)
    plot!(p5, t, [sol[network_idx[I_12_Q_IDX], i] for i in 1:length(t)], label="I_12_Q", linewith=2)
    plot!(p5, t, [sol[network_idx[I_23_D_IDX], i] for i in 1:length(t)], label="I_23_D", linewith=2)
    plot!(p5, t, [sol[network_idx[I_23_Q_IDX], i] for i in 1:length(t)], label="I_23_Q", linewith=2)
    plot!(p5, t, [sol[network_idx[I_13_D_IDX], i] for i in 1:length(t)], label="I_13_D", linewith=2)
    plot!(p5, t, [sol[network_idx[I_13_Q_IDX], i] for i in 1:length(t)], label="I_13_Q", linewith=2)
    savefig(p5, "../results/line_currrents.png")

    return sol
end

# Run the simulation
sol = run_simulation("../data/ThreeBusMultiLoad.raw")  # Adjust the path to your network file

println("Simulation complete!")