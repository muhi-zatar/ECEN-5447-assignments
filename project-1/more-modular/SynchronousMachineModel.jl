module SynchronousMachineModels

using DifferentialEquations
using Plots
using LinearAlgebra

export SynchronousMachine, AVR, TurbineGovernor, Network, PowerSystem
export simulate, analyze_results, apply_voltage_step, apply_load_change, apply_combined_perturbation

# =====================================================================
# Parameter Definitions
# =====================================================================

"""
    MachineParams

Parameters for the synchronous machine model.
"""
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
    Ra::Float64     # Armature resistance
end

"""
    AVRParams

Parameters for the IEEE Type-1 Automatic Voltage Regulator.
"""
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

"""
    TurbineGovParams

Parameters for the steam turbine and governor system.
"""
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

"""
    NetworkParams

Parameters for the external network connecting to the infinite bus.
"""
struct NetworkParams
    R_e::Float64    # External resistance
    X_e::Float64    # External reactance
    V_∞::Float64    # Infinite bus voltage
end

# =====================================================================
# Component Classes
# =====================================================================

"""
    SynchronousMachine

Represents a synchronous machine with its parameters and state.
"""
mutable struct SynchronousMachine
    params::MachineParams
    
    # State variables
    δ::Float64       # Rotor angle
    ω::Float64       # Rotor speed
    Eq_p::Float64    # q-axis transient voltage
    Ed_p::Float64    # d-axis transient voltage
    Eq_pp::Float64   # q-axis subtransient voltage
    Ed_pp::Float64   # d-axis subtransient voltage
    
    # Constructor with default parameters
    function SynchronousMachine()
        params = MachineParams(
            3.5,    # H
            2.0,    # D
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
            0.003   # Ra
        )
        
        # Initialize state variables to zero
        new(params, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    end
    
    # Constructor with custom parameters
    function SynchronousMachine(params::MachineParams)
        new(params, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    end
end

"""
    AVR

Represents an IEEE Type-1 Automatic Voltage Regulator.
"""
mutable struct AVR
    params::AVRParams
    
    # State variables
    Vr::Float64      # Regulator voltage
    Vf::Float64      # Field voltage
    Rf::Float64      # Rate feedback
    
    # Reference voltage
    Vref::Float64
    
    # Constructor with default parameters
    function AVR()
        params = AVRParams(
            200.0,  # Ka
            0.02,   # Ta
            1.0,    # Ke
            0.2,    # Te
            0.05,   # Kf
            1.0,    # Tf
            5.0,    # Vr_max
            -5.0    # Vr_min
        )
        
        # Initialize state variables to zero
        new(params, 0.0, 0.0, 0.0, 1.0)
    end
    
    # Constructor with custom parameters
    function AVR(params::AVRParams)
        new(params, 0.0, 0.0, 0.0, 1.0)
    end
end

"""
    TurbineGovernor

Represents a steam turbine with governor.
"""
mutable struct TurbineGovernor
    params::TurbineGovParams
    
    # State variables
    Pg::Float64      # Governor output
    Pm_hp::Float64   # High pressure turbine mechanical power
    Pm_lp::Float64   # Low pressure turbine mechanical power
    
    # Speed reference
    ω_ref::Float64
    
    # Constructor with default parameters
    function TurbineGovernor()
        params = TurbineGovParams(
            0.05,   # R
            0.2,    # Tg
            0.3,    # T_ch
            7.0,    # T_rh
            0.3,    # F_hp
            0.7,    # F_lp
            1.1,    # P_max
            0.0     # P_min
        )
        
        # Initialize state variables to zero
        new(params, 0.0, 0.0, 0.0, 1.0)
    end
    
    # Constructor with custom parameters
    function TurbineGovernor(params::TurbineGovParams)
        new(params, 0.0, 0.0, 0.0, 1.0)
    end
end

"""
    Network

Represents the external network connecting to the infinite bus.
"""
mutable struct Network
    params::NetworkParams
    
    # Constructor with default parameters
    function Network()
        params = NetworkParams(
            0.0,    # R_e
            0.5,    # X_e
            1.0     # V_∞
        )
        
        new(params)
    end
    
    # Constructor with custom parameters
    function Network(params::NetworkParams)
        new(params)
    end
end

"""
    PowerSystem

Complete power system with all components.
"""
mutable struct PowerSystem
    machine::SynchronousMachine
    avr::AVR
    turbine_gov::TurbineGovernor
    network::Network
    
    # Perturbation storage
    perturbations::Dict{Symbol, Tuple{Float64, Float64}}
    
    # Constructor with default components
    function PowerSystem()
        new(
            SynchronousMachine(),
            AVR(),
            TurbineGovernor(),
            Network(),
            Dict{Symbol, Tuple{Float64, Float64}}()
        )
    end
    
    # Constructor with custom components
    function PowerSystem(
        machine::SynchronousMachine,
        avr::AVR,
        turbine_gov::TurbineGovernor,
        network::Network
    )
        new(
            machine,
            avr,
            turbine_gov,
            network,
            Dict{Symbol, Tuple{Float64, Float64}}()
        )
    end
end

# =====================================================================
# State Management Functions
# =====================================================================

"""
    get_state_vector(system::PowerSystem)

Extract the complete state vector from the system components.
"""
function get_state_vector(system::PowerSystem)
    return [
        system.machine.δ,
        system.machine.ω,
        system.machine.Eq_p,
        system.machine.Ed_p,
        system.machine.Eq_pp,
        system.machine.Ed_pp,
        system.avr.Vr,
        system.avr.Vf,
        system.avr.Rf,
        system.turbine_gov.Pg,
        system.turbine_gov.Pm_hp,
        system.turbine_gov.Pm_lp
    ]
end

"""
    set_state_vector!(system::PowerSystem, u::Vector{Float64})

Update the system components from a state vector.
"""
function set_state_vector!(system::PowerSystem, u::Vector{Float64})
    system.machine.δ = u[1]
    system.machine.ω = u[2]
    system.machine.Eq_p = u[3]
    system.machine.Ed_p = u[4]
    system.machine.Eq_pp = u[5]
    system.machine.Ed_pp = u[6]
    system.avr.Vr = u[7]
    system.avr.Vf = u[8]
    system.avr.Rf = u[9]
    system.turbine_gov.Pg = u[10]
    system.turbine_gov.Pm_hp = u[11]
    system.turbine_gov.Pm_lp = u[12]
end

"""
    initialize!(system::PowerSystem; P0::Float64=0.8, Q0::Float64=0.0, Vt0::Float64=1.0)

Initialize the system to a steady state operating point.
"""
function initialize!(system::PowerSystem; P0::Float64=0.8, Q0::Float64=0.0, Vt0::Float64=1.0)
    # Extract parameters for convenience
    machine = system.machine
    avr = system.avr
    turbine_gov = system.turbine_gov
    network = system.network
    
    # Calculate initial rotor angle
    δ0 = asin(P0 * network.params.X_e / (Vt0 * network.params.V_∞))
    
    # Terminal voltages
    Vd0 = -network.params.V_∞ * sin(δ0)
    Vq0 = network.params.V_∞ * cos(δ0)
    
    # Calculate currents from desired power
    S0 = P0 + im*Q0
    I0 = conj(S0/Vt0)
    θ0 = angle(I0)
    
    Id0 = -abs(I0) * sin(θ0 + δ0)
    Iq0 = abs(I0) * cos(θ0 + δ0)
    
    # Calculate subtransient EMFs
    Eq_pp0 = Vd0 + machine.params.Ra * Id0 + machine.params.Xqpp * Iq0
    Ed_pp0 = Vq0 - machine.params.Ra * Iq0 - machine.params.Xdpp * Id0
    
    # Calculate transient EMFs
    Eq_p0 = Eq_pp0 + (machine.params.Xdp - machine.params.Xdpp) * Id0
    Ed_p0 = Ed_pp0 - (machine.params.Xqp - machine.params.Xqpp) * Iq0
    
    # Calculate required field voltage
    Vf0 = Eq_p0 + (machine.params.Xd - machine.params.Xdp) * Id0
    
    # AVR equilibrium
    Vr0 = Vf0 * avr.params.Ke
    Rf0 = avr.params.Kf * Vf0 / avr.params.Te
    
    # Governor equilibrium
    Pg0 = P0
    Pm_hp0 = turbine_gov.params.F_hp * P0
    Pm_lp0 = turbine_gov.params.F_lp * P0
    
    # Set initial state
    machine.δ = δ0
    machine.ω = 1.0
    machine.Eq_p = Eq_p0
    machine.Ed_p = Ed_p0
    machine.Eq_pp = Eq_pp0
    machine.Ed_pp = Ed_pp0
    
    avr.Vr = Vr0
    avr.Vf = Vf0
    avr.Rf = Rf0
    avr.Vref = 1.0
    
    turbine_gov.Pg = Pg0
    turbine_gov.Pm_hp = Pm_hp0
    turbine_gov.Pm_lp = Pm_lp0
    turbine_gov.ω_ref = 1.0
    
    # Validate initialization by checking derivatives
    u0 = get_state_vector(system)
    du0 = similar(u0)
    dynamics!(du0, u0, system, 0.0)
    
    # Print derivatives for validation
    println("Initial derivatives:")
    for i in 1:length(du0)
        println("du[$i] = $(du0[i])")
    end
    
    if maximum(abs.(du0)) > 1e-2
        @warn "System not at equilibrium. Maximum derivative: $(maximum(abs.(du0)))"
        
        # Optional: Run a short stabilization simulation with high damping
        # This could be implemented for better initialization
    end
    
    return system
end

# =====================================================================
# System Dynamics
# =====================================================================

"""
    dynamics!(du::Vector{Float64}, u::Vector{Float64}, system::PowerSystem, t::Float64)

Calculate the derivatives for all state variables.
"""
function dynamics!(du::Vector{Float64}, u::Vector{Float64}, system::PowerSystem, t::Float64)
    # Extract components for convenience
    machine = system.machine.params
    avr = system.avr.params
    turbine_gov = system.turbine_gov.params
    network = system.network.params
    perturbations = system.perturbations
    
    # Unpack state variables
    δ, ω, Eq_p, Ed_p, Eq_pp, Ed_pp, Vr, Vf, Rf, Pg, Pm_hp, Pm_lp = u
    
    # Calculate terminal voltages
    Vd = -network.V_∞ * sin(δ)
    Vq = network.V_∞ * cos(δ)
    
    # Matrix solution for currents
    A = [machine.Ra -machine.Xqpp; machine.Xdpp machine.Ra]
    b = [Eq_pp - Vd; Ed_pp + Vq]
    
    # Solve for currents
    currents = A \ b
    Id = currents[1]
    Iq = currents[2]
    
    # Calculate terminal voltage and power
    Vt = sqrt(Vd^2 + Vq^2)
    Pe = Vd*Id + Vq*Iq + machine.Ra*(Id^2 + Iq^2)
    Pm = Pm_hp + Pm_lp
    
    # Mechanical dynamics (swing equation)
    du[1] = ω - 1.0
    du[2] = (Pm - Pe - machine.D*(ω - 1.0))/(2.0*machine.H)
    
    # Electrical dynamics
    du[3] = (-Eq_p - (machine.Xd - machine.Xdp) * Id + Vf) / machine.Td0p
    du[4] = (-Ed_p + (machine.Xq - machine.Xqp) * Iq) / machine.Tq0p
    du[5] = (-Eq_pp + Eq_p + (machine.Xdp - machine.Xdpp) * Id) / machine.Td0pp
    du[6] = (-Ed_pp + Ed_p - (machine.Xqp - machine.Xqpp) * Iq) / machine.Tq0pp
    
    # Get reference voltage with any perturbations
    Vref = system.avr.Vref
    if haskey(perturbations, :v_ref)
        v_ref_change, v_ref_time = perturbations[:v_ref]
        if t >= v_ref_time
            Vref += v_ref_change
        end
    end
    
    # AVR dynamics
    du[7] = (-Vr + avr.Ka * (Vref - Vt - Rf)) / avr.Ta
    
    # Apply limits to regulator voltage
    if Vr > avr.Vr_max && du[7] > 0
        du[7] = 0.0
    elseif Vr < avr.Vr_min && du[7] < 0
        du[7] = 0.0
    end
    
    du[8] = (-Vf * avr.Ke + Vr) / avr.Te
    du[9] = (-Rf + avr.Kf * du[8]) / avr.Tf
    
    # Get reference speed with any perturbations
    ω_ref = system.turbine_gov.ω_ref
    if haskey(perturbations, :p_ref)
        p_ref_change, p_ref_time = perturbations[:p_ref]
        if t >= p_ref_time
            ω_ref += p_ref_change
        end
    end
    
    # Governor dynamics
    du[10] = (-Pg + (ω_ref - ω) / turbine_gov.R) / turbine_gov.Tg
    
    # Apply limits to governor output
    if Pg > turbine_gov.P_max && du[10] > 0
        du[10] = 0.0
    elseif Pg < turbine_gov.P_min && du[10] < 0
        du[10] = 0.0
    end
    
    # Turbine dynamics
    du[11] = (-Pm_hp + turbine_gov.F_hp * Pg) / turbine_gov.T_ch
    du[12] = (-Pm_lp + turbine_gov.F_lp * Pg + 
              (Pm_hp - turbine_gov.F_hp * Pg) / turbine_gov.T_rh) / turbine_gov.T_ch
    
    # Add numerical stability safeguards
    for i in 1:length(du)
        if abs(du[i]) > 100.0
            du[i] = sign(du[i]) * 100.0
        end
    end
end

# =====================================================================
# Simulation and Perturbation Functions
# =====================================================================

"""
    apply_voltage_step!(system::PowerSystem, step_size::Float64, step_time::Float64)

Apply a voltage reference step perturbation.
"""
function apply_voltage_step!(system::PowerSystem, step_size::Float64, step_time::Float64)
    system.perturbations[:v_ref] = (step_size, step_time)
    return system
end

"""
    apply_load_change!(system::PowerSystem, change_size::Float64, change_time::Float64)

Apply a load change perturbation.
"""
function apply_load_change!(system::PowerSystem, change_size::Float64, change_time::Float64)
    system.perturbations[:p_ref] = (change_size, change_time)
    return system
end

"""
    apply_combined_perturbation!(system::PowerSystem, v_step::Float64, v_time::Float64, 
                                p_change::Float64, p_time::Float64)

Apply both voltage step and load change perturbations.
"""
function apply_combined_perturbation!(system::PowerSystem, v_step::Float64, v_time::Float64, 
                                     p_change::Float64, p_time::Float64)
    system.perturbations[:v_ref] = (v_step, v_time)
    system.perturbations[:p_ref] = (p_change, p_time)
    return system
end

"""
    clear_perturbations!(system::PowerSystem)

Clear all perturbations from the system.
"""
function clear_perturbations!(system::PowerSystem)
    empty!(system.perturbations)
    return system
end

"""
    simulate(system::PowerSystem, tspan::Tuple{Float64, Float64})

Run a time-domain simulation of the power system.
"""
function simulate(system::PowerSystem, tspan::Tuple{Float64, Float64})
    # Get initial state vector
    u0 = get_state_vector(system)
    
    # Define ODE problem
    dynamics_wrapper = (du, u, p, t) -> dynamics!(du, u, p, t)
    prob = ODEProblem(dynamics_wrapper, u0, tspan, system)
    
    # Solve ODE problem
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
    
    return sol
end

# =====================================================================
# Analysis and Visualization
# =====================================================================

"""
    analyze_results(sol, system::PowerSystem, title_prefix::String="")

Analyze and plot simulation results.
"""
function analyze_results(sol, system::PowerSystem, title_prefix::String="")
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
    
    # Calculate terminal voltage and power at each time step
    Vt = zeros(length(sol.t))
    Pe = zeros(length(sol.t))
    
    for i in 1:length(sol.t)
        # Extract state for this time step
        state = sol.u[i]
        
        # Update system with this state
        set_state_vector!(system, state)
        
        # Calculate terminal variables
        Vd = -system.network.params.V_∞ * sin(state[1])
        Vq = system.network.params.V_∞ * cos(state[1])
        
        # Matrix solution for currents
        A = [system.machine.params.Ra -system.machine.params.Xqpp; 
             system.machine.params.Xdpp system.machine.params.Ra]
        b = [state[5] - Vd; state[6] + Vq]
        
        currents = A \ b
        Id = currents[1]
        Iq = currents[2]
        
        # Store results
        Vt[i] = sqrt(Vd^2 + Vq^2)
        Pe[i] = Vd*Id + Vq*Iq + system.machine.params.Ra*(Id^2 + Iq^2)
    end
    
    # Create plots
    p1 = plot(sol.t, ω, label="ω", title="$(title_prefix) Rotor Speed", 
              ylabel="ω (pu)", xlabel="Time (s)")
    
    p2 = plot(sol.t, δ * 180/π, label="δ", title="$(title_prefix) Rotor Angle", 
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
    combined_plot = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(800, 600), legend=true)
    
    return (combined_plot, (δ=δ, ω=ω, Vt=Vt, Vf=Vf, Pe=Pe, Pm=(Pm_hp, Pm_lp)))
end

# =====================================================================
# Convenience Functions for Perturbation Studies
# =====================================================================

"""
    voltage_reference_step(system::PowerSystem, step_size::Float64=0.05, step_time::Float64=1.0, 
                         sim_time::Float64=15.0)

Run a simulation with a voltage reference step perturbation.
"""
function voltage_reference_step(system::PowerSystem, step_size::Float64=0.05, step_time::Float64=1.0, 
                              sim_time::Float64=15.0)
    # Clear any existing perturbations
    clear_perturbations!(system)
    
    # Apply voltage step
    apply_voltage_step!(system, step_size, step_time)
    
    # Run simulation
    sol = simulate(system, (0.0, sim_time))
    
    # Analyze results
    return analyze_results(sol, system, "Voltage Reference Step")
end

"""
    load_change_step(system::PowerSystem, change_size::Float64=-0.1, change_time::Float64=1.0, 
                   sim_time::Float64=20.0)

Run a simulation with a load change perturbation.
"""
function load_change_step(system::PowerSystem, change_size::Float64=-0.1, change_time::Float64=1.0, 
                        sim_time::Float64=20.0)
    # Clear any existing perturbations
    clear_perturbations!(system)
    
    # Apply load change
    apply_load_change!(system, change_size, change_time)
    
    # Run simulation
    sol = simulate(system, (0.0, sim_time))
    
    # Analyze results
    return analyze_results(sol, system, "Load Change Step")
end

"""
    combined_perturbation(system::PowerSystem, v_step::Float64=0.03, v_time::Float64=1.0, 
                        p_change::Float64=-0.05, p_time::Float64=5.0, sim_time::Float64=25.0)

Run a simulation with both voltage step and load change perturbations.
"""
function combined_perturbation(system::PowerSystem, v_step::Float64=0.03, v_time::Float64=1.0, 
                             p_change::Float64=-0.05, p_time::Float64=5.0, sim_time::Float64=25.0)
    # Clear any existing perturbations
    clear_perturbations!(system)
    
    # Apply both perturbations
    apply_combined_perturbation!(system, v_step, v_time, p_change, p_time)
    
    # Run simulation
    sol = simulate(system, (0.0, sim_time))
    
    # Analyze results
    return analyze_results(sol, system, "Combined Perturbation")
end

end # module