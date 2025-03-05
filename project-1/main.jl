using DifferentialEquations
using Plots
using LinearAlgebra

"""
    Comprehensive Synchronous Machine Model based on Sauer-Pai approach
    Including:
    - 6th order synchronous machine model
    - IEEE Type-1 Automatic Voltage Regulator (AVR)
    - Steam Turbine with reheat
    - Governor with droop control
"""

# Machine Parameters (per unit)
struct MachineParams
    H::Float64      # Inertia constant (s)
    D::Float64      # Damping coefficient
    Xd::Float64     # d-axis synchronous reactance
    Xq::Float64     # q-axis synchronous reactance
    Xdp::Float64    # d-axis transient reactance
    Xqp::Float64    # q-axis transient reactance
    Xdpp::Float64   # d-axis subtransient reactance
    Xqpp::Float64   # q-axis subtransient reactance
    Xl::Float64     # Leakage reactance
    Td0p::Float64   # d-axis transient open-circuit time constant (s)
    Tq0p::Float64   # q-axis transient open-circuit time constant (s)
    Td0pp::Float64  # d-axis subtransient open-circuit time constant (s)
    Tq0pp::Float64  # q-axis subtransient open-circuit time constant (s)
    Taa::Float64    # d-axis additional leakage time constant (s)
    Ra::Float64     # Armature resistance
end

# IEEE Type-1 AVR Parameters
struct AVRParams
    Ka::Float64     # Amplifier gain
    Ta::Float64     # Amplifier time constant (s)
    Ke::Float64     # Exciter constant
    Te::Float64     # Exciter time constant (s)
    Kf::Float64     # Stabilizer gain
    Tf::Float64     # Stabilizer time constant (s)
    Vr_max::Float64 # Maximum regulator voltage
    Vr_min::Float64 # Minimum regulator voltage
end

# Turbine-Governor Parameters
struct TurbineGovParams
    R::Float64      # Droop coefficient
    Tg::Float64     # Governor time constant (s)
    T_ch::Float64   # Steam chest time constant (s)
    T_rh::Float64   # Reheat time constant (s)
    F_hp::Float64   # High pressure turbine fraction
    F_lp::Float64   # Low pressure turbine fraction
    P_max::Float64  # Maximum power
    P_min::Float64  # Minimum power
end

# Network Parameters
struct NetworkParams
    R_e::Float64    # External resistance
    X_e::Float64    # External reactance
    V_∞::Float64    # Infinite bus voltage
end

# Define default parameters for a typical machine
function default_machine_params()
    return MachineParams(
        2.0,    # H
        10.0,    # D
        1.81,   # Xd
        1.76,   # Xq
        0.3,    # Xdp
        0.65,   # Xqp
        0.23,   # Xdpp
        0.25,   # Xqpp
        0.15,   # Xl
        8.0,    # Td0p
        0.4,    # Tq0p
        0.03,   # Td0pp
        0.05,   # Tq0pp
        0.00,   # Taa
        0.003   # Ra
    )
end

# Define default parameters for IEEE Type-1 AVR
function default_avr_params()
    return AVRParams(
        200.0,  # Ka
        0.02,   # Ta
        1.0,    # Ke
        0.2,    # Te
        0.05,   # Kf
        1.0,    # Tf
        5.0,    # Vr_max
        -5.0    # Vr_min
    )
end

# Define default parameters for turbine and governor
function default_turbine_gov_params()
    return TurbineGovParams(
        0.05,   # R
        0.2,    # Tg
        0.3,    # T_ch
        7.0,    # T_rh
        0.3,    # F_hp
        0.7,    # F_lp
        1.1,    # P_max
        0.0     # P_min
    )
end

# Define default network parameters
function default_network_params()
    return NetworkParams(
        0.0,    # R_e
        0.5,    # X_e
        1.0     # V_∞
    )
end

"""
    Synchronous machine model with AVR, turbine and governor
    State variables:
    u[1] = δ       : rotor angle
    u[2] = ω       : rotor speed
    u[3] = E'q     : q-axis transient voltage
    u[4] = E'd     : d-axis transient voltage
    u[5] = E''q    : q-axis subtransient voltage
    u[6] = E''d    : d-axis subtransient voltage
    u[7] = Vr      : regulator voltage
    u[8] = Vf      : exciter output
    u[9] = Rf      : rate feedback
    u[10] = Pg     : governor output
    u[11] = Pm_hp  : High pressure turbine mechanical power
    u[12] = Pm_lp  : Low pressure turbine mechanical power
"""
function synchronous_machine_dynamics!(du, u, p, t)
    machine, avr, turbine_gov, network, perturbations = p

    # Unpack state variables
    δ, ω, Eq_p, Ed_p, ψq_pp, ψd_pp, Vr, Vf, Rf, Pg, Pm_hp, Pm_lp = u

    # Calculate terminal conditions
    # Stator voltage equations
    Vd = network.V_∞ * sin(δ)       # Eqn 15.4
    Vq = network.V_∞ * cos(δ)       # Eqn 15.4

    # Terminal voltage
    Vt = sqrt(Vd^2 + Vq^2)

    # Calculate machine currents
    # Helper variables (Eqn 15.14)
    γd1 = (machine.Xdpp - machine.Xl) / (machine.Xdp - machine.Xl)
    γq1 = (machine.Xqpp - machine.Xl) / (machine.Xqp - machine.Xl)
    γd2 = (1 - γd1) / (machine.Xdp - machine.Xl)
    γq2 = (1 - γq1) / (machine.Xqp - machine.Xl)

    # Current equations (Use 15.11 and 15.15 to eliminate ψq and ψd, solve for Id and Iq)
    A = [machine.Xdpp machine.Ra; machine.Ra -machine.Xqpp]
    b = [-Vq + γd1 * Eq_p + (1 - γd1) * ψd_pp; -Vd + γq1 * Ed_p + (1 - γq1) * ψq_pp]
    currents = A \ b
    Id = currents[1]
    Iq = currents[2]

    # Calculate shaft conditions
    # Electromagnetic torque/power (Use 15.11 and 15.6 to eliminate ψq and ψd, solve for Te)
    Pe = Vd * Id + Vq * Iq + machine.Ra * (Id^2 + Iq^2)
    # Te = Pe   # in per unit, torque equals power

    # Total mechanical power input
    Pm = Pm_hp + Pm_lp

    # Synchronous speed in pu
    ω_s = 1.0

    # Shaft mechanical equations (Eqn 15.5)
    du[1] = ω - ω_s  # dδ/dt
    du[2] = (Pm - Pe - machine.D * (ω - ω_s)) / (2.0 * machine.H)

    # Electrical dynamics
    du[3] = (-Eq_p - (machine.Xd - machine.Xdp) * (Id - γd2 * ψd_pp - (1 - γd1) * Id + γd2 * Eq_p) + Vf) / machine.Td0p    # dE'q/dt
    du[4] = (-Ed_p + (machine.Xq - machine.Xqp) * (Iq - γq2 * ψq_pp - (1 - γq1) * Iq - γd2 * Ed_p)) / machine.Tq0p         # dE'd/dt
    du[5] = (-ψd_pp + Eq_p - (machine.Xdp - machine.Xl) * Id) / machine.Td0pp                                              # dE''q/dt
    du[6] = (-ψq_pp - Ed_p - (machine.Xqp - machine.Xl) * Iq) / machine.Tq0pp                                              # dE''d/dt

    # Apply voltage perturbation if specified
    Vref = 1.0
    if haskey(perturbations, :v_ref)
        v_ref_change, v_ref_time = perturbations[:v_ref]
        if t >= v_ref_time
            Vref += v_ref_change
        end
    end

    # AVR dynamics (IEEE Type-1)
    du[7] = (-Vr + avr.Ka * (Vref - Vt - Rf)) / avr.Ta  # dVr/dt
    # Limit Vr within bounds
    if Vr > avr.Vr_max && du[7] > 0
        du[7] = 0.0
    elseif Vr < avr.Vr_min && du[7] < 0
        du[7] = 0.0
    end

    du[8] = (-Vf * avr.Ke + Vr) / avr.Te  # dVf/dt
    du[9] = (-Rf + avr.Kf * (du[8])) / avr.Tf  # dRf/dt

    # Apply load perturbation if specified
    ω_ref = 1.0
    # Apply frequency reference perturbation if defined
    if haskey(perturbations, :p_ref)
        p_ref_change, p_ref_time = perturbations[:p_ref]
        if t >= p_ref_time
            ω_ref += p_ref_change
        end
    end

    # Governor dynamics
    du[10] = (-Pg + (ω_ref - ω) / turbine_gov.R) / turbine_gov.Tg           # dPg/dt
    # Limit Pg within bounds
    if Pg > turbine_gov.P_max && du[10] > 0
        du[10] = 0.0
    elseif Pg < turbine_gov.P_min && du[10] < 0
        du[10] = 0.0
    end

    # Turbine dynamics
    du[11] = (-Pm_hp + turbine_gov.F_hp * Pg) / turbine_gov.T_ch            # dPm_hp/dt
    du[12] = (-Pm_lp + turbine_gov.F_lp * Pg + (Pm_hp - turbine_gov.F_hp * Pg) / turbine_gov.T_rh) / turbine_gov.T_ch  # dPm_lp/dt

    # Add numerical stability safeguards
    for i in 1:length(du)
        # Limit the rate of change to reasonable values
        if abs(du[i]) > 100.0
            du[i] = sign(du[i]) * 100.0
        end
    end
end

#=

# Function to set up the initial conditions for steady state operation
function init_steady_state(machine, avr, turbine_gov, network)
    # Assume initial power output P0 and terminal voltage Vt0
    P0 = 0.8
    Vt0 = 1.0

    # Calculate initial machine angle for this power output
    δ0 = asin(P0 * network.X_e / (Vt0 * network.V_∞))

    # Calculate initial Vd and Vq
    Vd0 = -network.V_∞ * sin(δ0)
    Vq0 = network.V_∞ * cos(δ0)

    S0 = P0  # Assuming unity power factor for simplicity
    I0 = S0 / Vt0
    θ0 = 0.0  # Power angle (assuming unity power factor)
    # Calculate initial Id and Iq
    Id0 = -I0 * sin(θ0 + δ0)
    Iq0 = I0 * cos(θ0 + δ0)
    # Id0 = (P0 * Vd0) / (Vd0^2 + Vq0^2) # Caused inequilibiruim 
    # Iq0 = (P0 * Vq0) / (Vd0^2 + Vq0^2) # Caused inequilibiruim

    # Calculate initial internal voltages
    Eq_pp0 = Vd0 + machine.Ra * Id0 + machine.Xdpp * Iq0
    Ed_pp0 = Vq0 - machine.Ra * Iq0 - machine.Xqpp * Id0

    # Calculate intermediate variables
    Ψ1d0 = Eq_pp0 + (machine.Xdpp - machine.Xl) * Id0
    Ψ1q0 = Ed_pp0 - (machine.Xqpp - machine.Xl) * Iq0

    # Calculate transient voltages
    Eq_p0 = Eq_pp0 + (machine.Xdp - machine.Xdpp) * Id0
    Ed_p0 = Ed_pp0 - (machine.Xqp - machine.Xqpp) * Iq0

    # Initial exciter output
    Vf0 = Eq_p0 + (machine.Xd - machine.Xdp) * Id0

    # Initial AVR signals
    Rf0 = 0.0  # assume steady state
    Vr0 = Vf0 * avr.Ke  # steady state value that account for exciter constant

    # Initial governor and turbine outputs
    Pg0 = P0
    Pm_hp0 = turbine_gov.F_hp * P0
    Pm_lp0 = turbine_gov.F_lp * P0

    u0 = [δ0, 1.0, Eq_p0, Ed_p0, Eq_pp0, Ed_pp0, Vr0, Vf0, Rf0, Pg0, Pm_hp0, Pm_lp0]

    # The following lines are sanity checks for the initial conditions
    machine_high_damping = MachineParams(
        machine.H, 10.0, machine.Xd, machine.Xq, machine.Xdp, machine.Xqp,
        machine.Xdpp, machine.Xqpp, machine.Xl, machine.Td0p, machine.Tq0p,
        machine.Td0pp, machine.Tq0pp, machine.Ra
    )

    p_init = (machine_high_damping, avr, turbine_gov, network, Dict())
    prob_init = ODEProblem(synchronous_machine_dynamics!, u0, (0.0, 5.0), p_init)
    sol_init = solve(prob_init, Tsit5(), reltol=1e-6, abstol=1e-6)
    # Done checking

    # Return initial state vector
    return sol_init.u[end]
end
=#

function init_steady_state(machine, avr, turbine_gov, network)
    # Start with a simpler power flow solution
    P0 = 0.8  # Active power output
    Q0 = 0.0  # Reactive power (assuming unity power factor)
    Vt0 = 1.0  # Terminal voltage

    # Network calculations
    V_inf = network.V_∞
    R_e = network.R_e
    X_e = network.X_e
    Z_e = sqrt(R_e^2 + X_e^2)
    theta_e = atan(X_e, R_e)

    δ0 = asin(P0 * network.X_e / network.V_∞)

    # Calculate machine angle using power flow equations
    # changed name
    # delta0 = asin((P0*Z_e)/(Vt0*V_inf)) + theta_e

    # Terminal voltages
    # Vd0 = -V_inf*sin(delta0)
    # Vq0 = V_inf*cos(delta0)
    Vd0 = -network.V_∞ * sin(δ0)
    Vq0 = network.V_∞ * cos(δ0)

    # Machine currents (simplified network equations)
    S0 = complex(P0, Q0)
    I0_mag = abs(S0 / Vt0)
    I0_ang = angle(S0 / Vt0) # Power factor angle

    Id0 = (P0 * Vd0 + Q0 * Vq0) / (Vd0^2 + Vq0^2)
    Iq0 = (P0 * Vq0 - Q0 * Vd0) / (Vd0^2 + Vq0^2)

    # Start with simplified model steady state
    # Calculate subtransient EMFs
    Eq_pp0 = Vd0 + machine.Ra * Id0 + machine.Xqpp * Iq0
    Ed_pp0 = Vq0 - machine.Ra * Iq0 - machine.Xdpp * Id0

    # Calculate transient EMFs
    Eq_p0 = Eq_pp0 + (machine.Xdp - machine.Xdpp) * Id0
    Ed_p0 = Ed_pp0 - (machine.Xqp - machine.Xqpp) * Iq0

    # Field voltage needed for this operating point
    Vf0 = Eq_p0 + (machine.Xd - machine.Xdp) * Id0

    # Initial governor and exciter values
    Vr0 = Vf0 * avr.Ke # AVR equilibrium - set Vr to match Vf

    Rf0 = avr.Kf * Vf0 / avr.Te #  Rate feedback in equilibrium
    Pg0 = P0
    Pm_hp0 = turbine_gov.F_hp * P0
    Pm_lp0 = turbine_gov.F_lp * P0

    # Return true steady-state initial conditions
    return [δ0, 1.0, Eq_p0, Ed_p0, Eq_pp0, Ed_pp0, Vr0, Vf0, Rf0, Pg0, Pm_hp0, Pm_lp0]
end

# Function to run a simulation with specified perturbations
function run_simulation(perturbations, tspan=(0.0, 10.0))
    # Set up parameters
    machine = default_machine_params()
    avr = default_avr_params()
    turbine_gov = default_turbine_gov_params()
    network = default_network_params()

    # Pack parameters for ODE solver
    p = (machine=machine, avr=avr, turbine_gov=turbine_gov, network=network, perturbations=perturbations)

    function f(du, u, p, t)::Nothing
        synchronous_machine_dynamics!(du, u, (p.machine, p.avr, p.turbine_gov, p.network, p.perturbations), t)
        nothing  # Explicit return of nothing (this was to solve an error, not sure I understand it)
    end
    # Set up initial conditions
    u0 = init_steady_state(machine, avr, turbine_gov, network)

    # Define ODE problem
    prob = ODEProblem(f, u0, tspan, p)

    # Solve ODE problem
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

    # INITIAL PARAMETERS SANITY CHECK
    u0 = init_steady_state(machine, avr, turbine_gov, network)
    du0 = zeros(12)
    synchronous_machine_dynamics!(du0, u0, (machine, avr, turbine_gov, network, Dict()), 0.0)

    println("Initial derivatives:")
    for i in 1:length(du0)
        println("du[$i] = $(du0[i])")
    end

    if maximum(abs.(du0)) > 1e-3
        @warn "System not at equilibrium. Maximum derivative: $(maximum(abs.(du0)))"
    end

    return sol
end

# Function to analyze and plot results
function analyze_results(sol, title_prefix="")
    # Extract state variables
    δ = [u[1] for u in sol.u]
    ω = [u[2] for u in sol.u]
    Eq_p = [u[3] for u in sol.u]
    Ed_p = [u[4] for u in sol.u]
    Eq_pp = [u[5] for u in sol.u]
    Ed_pp = [u[6] for u in sol.u]
    Vr = [u[7] for u in sol.u]
    Vf = [u[8] for u in sol.u]
    Pg = [u[10] for u in sol.u]
    Pm_hp = [u[11] for u in sol.u]
    Pm_lp = [u[12] for u in sol.u]

    # Calculate terminal voltage at each time step
    Vt = similar(δ)
    Pe = similar(δ)

    machine = default_machine_params()
    network = default_network_params()

    for i in 1:length(sol.t)
        # Terminal conditions
        Vd = -network.V_∞ * sin(δ[i])
        Vq = network.V_∞ * cos(δ[i])

        # Calculate terminal voltage
        Vt[i] = sqrt(Vd^2 + Vq^2)

        # Calculate currents
        Id = (Eq_pp[i] - Vd - machine.Ra * (Vq - Ed_pp[i]) / machine.Xqpp) /
             (machine.Ra^2 / machine.Xqpp + machine.Xdpp)
        Iq = (Ed_pp[i] + Vq - machine.Ra * (Eq_pp[i] - Vd) / machine.Xdpp) /
             (machine.Ra^2 / machine.Xdpp + machine.Xqpp)

        # Calculate electrical power
        Pe[i] = (Vd + machine.Ra * Id) * Id + (Vq + machine.Ra * Iq) * Iq
    end

    # Plot results
    p1 = plot(sol.t, ω, label="ω", title="$(title_prefix) Rotor Speed",
        ylabel="ω (pu)", xlabel="Time (s)")

    p2 = plot(sol.t, δ * 180 / π, label="δ", title="$(title_prefix) Rotor Angle",
        ylabel="δ (degrees)", xlabel="Time (s)")

    p3 = plot(sol.t, Vt, label="Vt", title="$(title_prefix) Terminal Voltage",
        ylabel="Voltage (pu)", xlabel="Time (s)")

    p4 = plot(sol.t, Vf, label="Vf", title="$(title_prefix) Field Voltage",
        ylabel="Vf (pu)", xlabel="Time (s)")

    p5 = plot(sol.t, Pe, label="Pe", title="$(title_prefix) Electrical Power",
        ylabel="Power (pu)", xlabel="Time (s)")

    p6 = plot(sol.t, [Pm_hp Pm_lp], label=["Pm_hp" "Pm_lp"],
        title="$(title_prefix) Mechanical Power",
        ylabel="Power (pu)", xlabel="Time (s)")

    # Combine plots
    plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(800, 600), legend=true)
end

# Run simulations with different perturbations

# 1. Step change in reference voltage (AVR test)
function voltage_reference_step()
    perturbations = Dict(:v_ref => (0.05, 1.0))  # 5% increase at t=1s
    sol = run_simulation(perturbations, (0.0, 15.0))
    analyze_results(sol, "Voltage Reference Step")
    return sol
end

# 2. Load change (Governor test)
function load_change_step()
    perturbations = Dict(:p_ref => (-0.1, 1.0))  # 10% load increase at t=1s
    sol = run_simulation(perturbations, (0.0, 20.0))
    analyze_results(sol, "Load Change Step")
    return sol
end

# 3. Combined perturbation
function combined_perturbation()
    perturbations = Dict(
        :v_ref => (0.03, 1.0),   # 3% voltage increase at t=1s
        :p_ref => (-0.05, 5.0)   # 5% load increase at t=5s
    )
    sol = run_simulation(perturbations, (0.0, 25.0))
    analyze_results(sol, "Combined Perturbation")
    return sol
end

# Main function to run all simulations
function main()
    println("Running Voltage Reference Step Test...")
    sol1 = voltage_reference_step()
    savefig("voltage_step_response.png")

    println("Running Load Change Test...")
    sol2 = load_change_step()
    savefig("load_change_response.png")

    println("Running Combined Perturbation Test...")
    sol3 = combined_perturbation()
    savefig("combined_perturbation_response.png")

    return sol1, sol2, sol3
end

# Execute the main function
main()