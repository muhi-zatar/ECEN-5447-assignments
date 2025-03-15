# Running simulation for generator and shaft in steady state.

using Plots
using DifferentialEquations
using LinearAlgebra

# Include models
include("SauerPaiMachineModel.jl")
include("PowerFlowWrapper.jl")

using .SauerPaiMachineModel
using .PowerFlowWrapper

"""
Simulates a synchronous machine with the Sauer-Pai model in steady state.
"""
function run_simulation(network_file_path, simulation_time=5.0, dt=0.01)

    v_mag, v_angle, P, Q = get_powerflow_results(network_file_path)
    V_terminal = v_mag * exp(im * v_angle)

    println("Power flow results:")
    println("  Terminal voltage magnitude: $v_mag pu")
    println("  Terminal voltage angle: $(v_angle * 180/π) degrees")
    println("  Active power: $P pu")
    println("  Reactive power: $Q pu")

    # Create machine with given parameters
    machine = SauerPaiMachine(
        R=0.002,
        X_d=1.79,
        X_q=1.71,
        Xd_p=0.169,
        Xq_p=0.228,
        Xd_pp=0.135,
        Xq_pp=0.2,
        Xl=0.13,
        Td0_p=4.3,
        Tq0_p=0.85,
        Td0_pp=0.032,
        Tq0_pp=0.05,
        H=3.148,
        D=2.0,
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

    initial_states, V_f_init, τ_m_init = initialize_machine(machine, V_terminal, v_angle, P, Q)

    println("Initial machine states:")
    state_names = ["Delta", "Omega", "Eq_p", "Ed_p", "Psi_d_pp", "Psi_q_pp", "Torque_e", "Vf"]
    for (i, name) in enumerate(state_names)
        println("  $name: $(initial_states[i])")
    end

    function ode_system!(du, u, p, t)
        current_v_mag = v_mag
        current_v_angle = v_angle

        current_V_terminal = current_v_mag * exp(im * current_v_angle)

        update_machine_states!(u, du, current_V_terminal, machine)
    end

    tspan = (0.0, simulation_time)
    prob = ODEProblem(ode_system!, initial_states, tspan)

    sol = solve(prob, Tsit5(), saveat=dt)

    t = sol.t
    delta = [sol[i][1] for i in 1:length(t)]
    omega = [sol[i][2] for i in 1:length(t)]
    eq_p = [sol[i][3] for i in 1:length(t)]
    ed_p = [sol[i][4] for i in 1:length(t)]
    psi_d_pp = [sol[i][5] for i in 1:length(t)]
    psi_q_pp = [sol[i][6] for i in 1:length(t)]
    torque_e = [sol[i][7] for i in 1:length(t)]
    vf = [sol[i][8] for i in 1:length(t)]

    v_terminal_mag = ones(length(t)) .* v_mag
    v_terminal_angle = ones(length(t)) .* v_angle

    return (
        time=t,
        delta=delta,
        omega=omega,
        eq_p=eq_p,
        ed_p=ed_p,
        psi_d_pp=psi_d_pp,
        psi_q_pp=psi_q_pp,
        torque_e=torque_e,
        vf=vf,
        v_terminal_mag=v_terminal_mag,
        v_terminal_angle=v_terminal_angle,
        solution=sol
    )
end

function plot_results(results)

    # Plot boltage
    p1 = plot(results.time, results.v_terminal_mag,
        label="Terminal Voltage (pu)",
        title="Terminal Voltage Magnitude",
        xlabel="Time (s)",
        ylabel="Voltage (pu)",
        linewidth=2,
        grid=true,
        legend=:bottomright)

    savefig(p1, "voltage_ss_gen_shaft.png")

    # Plot rotor angle
    p2 = plot(results.time, results.delta .* (180 / π),
        label="Rotor Angle (deg)",
        title="Rotor Angle",
        xlabel="Time (s)",
        ylabel="Angle (degrees)",
        linewidth=2,
        grid=true,
        legend=:bottomright,
        ylims=(results.delta[1] * (180 / π) - 0.1, results.delta[1] * (180 / π) + 0.1)) # I wanted to see the changes more clearly on the graph (the oscillations).

    savefig(p2, "rotor_angle_ss_gen_shaft.png")

    # Plot rotor speed
    p3 = plot(results.time, results.omega,
        label="Rotor Speed (pu)",
        title="Rotor Speed",
        xlabel="Time (s)",
        ylabel="Speed (pu)",
        linewidth=2,
        grid=true,
        legend=:bottomright,
        ylims=(0.999, 1.001))

    savefig(p3, "rotor_speed_ss_gen_shaft.png")

    return p1, p2, p3
end

"""
Main function to run the simulation and plot results
"""
function main()
    network_file = "../data/ThreeBusMultiLoad.raw"

    # Run simulation
    results = run_simulation(network_file)

    p1, p2, p3 = plot_results(results)

    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end