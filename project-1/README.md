# Synchronous Machine Model with AVR, Turbine, and Governor

This README provides a detailed explanation of the synchronous machine model implemented in Julia using the DifferentialEquations package. The model follows the Sauer-Pai approach for power system dynamics and includes:

- 6th order synchronous machine model
- IEEE Type-1 Automatic Voltage Regulator (AVR)
- Steam turbine with reheat
- Governor with droop control

## Table of Contents
- [Required Packages](#required-packages)
- [Parameter Definitions](#parameter-definitions)
- [Dynamic Model](#dynamic-model)
- [Initialization](#initialization)
- [Running Simulations](#running-simulations)
- [Analyzing Results](#analyzing-results)

## Required Packages

```julia
using DifferentialEquations
using Plots
using LinearAlgebra
```

- `DifferentialEquations`: Provides solvers for the differential equations representing the machine dynamics
- `Plots`: Used for visualizing simulation results
- `LinearAlgebra`: Required for matrix operations in current calculations

## Parameter Definitions

### Machine Parameters

```julia
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
```

This struct defines the synchronous machine parameters:
- `H`: Inertia constant (in seconds) representing the stored kinetic energy at rated speed
- `D`: Damping coefficient representing mechanical damping torque proportional to speed deviation
- `Xd`, `Xq`: Synchronous reactances in d-axis and q-axis
- `Xdp`, `Xqp`: Transient reactances in d-axis and q-axis
- `Xdpp`, `Xqpp`: Subtransient reactances in d-axis and q-axis
- `Xl`: Leakage reactance of the stator winding
- `Td0p`, `Tq0p`: Open-circuit transient time constants
- `Td0pp`, `Tq0pp`: Open-circuit subtransient time constants
- `Ra`: Armature resistance

### AVR Parameters

```julia
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
```

This defines the IEEE Type-1 AVR parameters:
- `Ka`: Amplifier gain
- `Ta`: Amplifier time constant
- `Ke`: Exciter gain constant
- `Te`: Exciter time constant
- `Kf`: Stabilizer gain
- `Tf`: Stabilizer time constant
- `Vr_max`, `Vr_min`: Voltage regulator output limits

### Turbine-Governor Parameters

```julia
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
```

This defines turbine and governor parameters:
- `R`: Droop coefficient (pu frequency change per pu power change)
- `Tg`: Governor time constant
- `T_ch`: Steam chest time constant
- `T_rh`: Reheat time constant
- `F_hp`, `F_lp`: Power fractions for high and low pressure turbine sections
- `P_max`, `P_min`: Governor output power limits

### Network Parameters

```julia
struct NetworkParams
    R_e::Float64    # External resistance
    X_e::Float64    # External reactance
    V_∞::Float64    # Infinite bus voltage
end
```

This defines the external network parameters:
- `R_e`: External resistance connecting the machine to the infinite bus
- `X_e`: External reactance connecting the machine to the infinite bus
- `V_∞`: Infinite bus voltage magnitude

### Default Parameters

```julia
function default_machine_params()
    return MachineParams(
        3.5,    # H
        2.0,    # D (note: increased from 0.0 for stability)
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
end
```

This function creates default machine parameters for a typical generator. Note that damping (D) is set to 2.0 for numerical stability. In reality, it might be smaller.

Similar default parameter functions exist for the AVR, turbine-governor, and network.

## Dynamic Model

The synchronous machine model is defined with 12 state variables:

```julia
function synchronous_machine_dynamics!(du, u, p, t)
    # State variables:
    # u[1] = δ       : rotor angle
    # u[2] = ω       : rotor speed
    # u[3] = E'q     : q-axis transient voltage
    # u[4] = E'd     : d-axis transient voltage
    # u[5] = E''q    : q-axis subtransient voltage
    # u[6] = E''d    : d-axis subtransient voltage
    # u[7] = Vr      : regulator voltage
    # u[8] = Vf      : exciter output
    # u[9] = Rf      : rate feedback
    # u[10] = Pg     : governor output
    # u[11] = Pm_hp  : High pressure turbine mechanical power
    # u[12] = Pm_lp  : Low pressure turbine mechanical power
    
    machine, avr, turbine_gov, network, perturbations = p
    
    # Unpack state variables
    δ, ω, Eq_p, Ed_p, Eq_pp, Ed_pp, Vr, Vf, Rf, Pg, Pm_hp, Pm_lp = u
```

Here, we define the function that calculates derivatives (du) for all state variables. The parameters `p` are unpacked to access machine, AVR, turbine-governor, and network parameters, along with any perturbations.

### Network Interface

```julia
    # Calculate terminal voltages relative to infinite bus
    Vd = -network.V_∞ * sin(δ)
    Vq = network.V_∞ * cos(δ)
```

These lines transform the infinite bus voltage to the machine's rotor reference frame. The machine's rotor angle δ defines the angle between the q-axis and the infinite bus voltage.

### Current Calculation

```julia
    # Matrix solution for currents
    A = [machine.Ra -machine.Xqpp; machine.Xdpp machine.Ra]
    b = [Eq_pp - Vd; Ed_pp + Vq]
    currents = A \ b
    Id = currents[1]
    Iq = currents[2]
```

This solves the algebraic loop for stator currents using matrix operations. The stator voltage equations are:
- Vd + Ra*Id - Eq_pp + Xq_pp*Iq = 0
- Vq + Ra*Iq - Ed_pp - Xd_pp*Id = 0

We solve this 2x2 system to find Id and Iq.

### Power Calculation

```julia
    # Calculate terminal voltage magnitude
    Vt = sqrt(Vd^2 + Vq^2)
    
    # Calculate electrical power
    Pe = Vd*Id + Vq*Iq + machine.Ra*(Id^2 + Iq^2)
    
    # Mechanical power
    Pm = Pm_hp + Pm_lp
```

This calculates the terminal voltage magnitude and electrical power output. The mechanical power is the sum of high-pressure and low-pressure turbine outputs.

### Mechanical Dynamics (Swing Equation)

```julia
    # Mechanical dynamics
    ω_s = 1.0  # Synchronous speed in pu
    du[1] = ω - ω_s  # dδ/dt
    du[2] = (Pm - Pe - machine.D*(ω - ω_s))/(2.0*machine.H)  # dω/dt
```

These are the swing equations:
- dδ/dt = ω - ω_s (where ω_s = 1.0 pu)
- dω/dt = (Pm - Pe - D*Δω)/(2H)

The first equation relates rotor angle change to speed deviation. The second is Newton's law for rotational motion, where the torque accelerating the rotor is proportional to power imbalance.

### Electrical Dynamics

```julia
    # Electrical dynamics
    du[3] = (-Eq_p - (machine.Xd - machine.Xdp) * Id + Vf) / machine.Td0p    # dE'q/dt
    du[4] = (-Ed_p + (machine.Xq - machine.Xqp) * Iq) / machine.Tq0p         # dE'd/dt
    du[5] = (-Eq_pp + Eq_p + (machine.Xdp - machine.Xdpp) * Id) / machine.Td0pp  # dE''q/dt
    du[6] = (-Ed_pp + Ed_p - (machine.Xqp - machine.Xqpp) * Iq) / machine.Tq0pp  # dE''d/dt
```

These equations model the dynamics of the internal voltage variables:
- dE'q/dt: q-axis transient voltage, affected by field voltage Vf
- dE'd/dt: d-axis transient voltage
- dE''q/dt: q-axis subtransient voltage
- dE''d/dt: d-axis subtransient voltage

These represent the electromagnetic dynamics in the rotor circuits.

### AVR Dynamics

```julia
    # Apply reference voltage perturbation if defined
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
    
    # Exciter equation
    du[8] = (-Vf * avr.Ke + Vr) / avr.Te  # dVf/dt
    
    # Rate feedback
    du[9] = (-Rf + avr.Kf * (du[8])) / avr.Tf  # dRf/dt
```

This implements the IEEE Type-1 AVR model:
1. First checks for any voltage reference perturbations
2. Calculates regulator voltage dynamics with amplifier gain and time constant
3. Applies limits to regulator voltage
4. Models exciter dynamics that convert regulator voltage to field voltage
5. Implements rate feedback for stabilization

### Governor and Turbine Dynamics

```julia
    # Apply frequency reference perturbation if defined
    ω_ref = 1.0
    if haskey(perturbations, :p_ref)
        p_ref_change, p_ref_time = perturbations[:p_ref]
        if t >= p_ref_time
            ω_ref += p_ref_change
        end
    end
    
    # Governor dynamics
    du[10] = (-Pg + (ω_ref - ω) / turbine_gov.R) / turbine_gov.Tg  # dPg/dt
    
    # Limit Pg within bounds
    if Pg > turbine_gov.P_max && du[10] > 0
        du[10] = 0.0
    elseif Pg < turbine_gov.P_min && du[10] < 0
        du[10] = 0.0
    end
    
    # Turbine dynamics
    du[11] = (-Pm_hp + turbine_gov.F_hp * Pg) / turbine_gov.T_ch  # dPm_hp/dt
    du[12] = (-Pm_lp + turbine_gov.F_lp * Pg + (Pm_hp - turbine_gov.F_hp * Pg) / turbine_gov.T_rh) / turbine_gov.T_ch  # dPm_lp/dt
```

This models:
1. Governor response to frequency deviations (with droop R)
2. Power output limits
3. High-pressure turbine dynamics
4. Low-pressure turbine dynamics with reheat

The low-pressure turbine receives both direct power from the governor and reheated steam from the high-pressure turbine.

### Numerical Stability Safeguards

```julia
    # Add numerical stability safeguards
    for i in 1:length(du)
        # Limit the rate of change to reasonable values
        if abs(du[i]) > 100.0
            du[i] = sign(du[i]) * 100.0
        end
    end
```

This limits extreme derivative values to prevent numerical instability during integration.

## Initialization

Proper initialization is critical for a stable simulation. The `init_steady_state` function calculates initial conditions that satisfy the system's equilibrium equations.

```julia
function init_steady_state(machine, avr, turbine_gov, network)
    # Start with desired operating point
    P0 = 0.8
    Q0 = 0.0
    Vt0 = 1.0
    
    # Calculate initial angle
    δ0 = asin(P0 * network.X_e / network.V_∞)
    
    # Terminal voltages
    Vd0 = -network.V_∞ * sin(δ0)
    Vq0 = network.V_∞ * cos(δ0)
```

This begins by setting the desired operating point (active power, reactive power, and terminal voltage) and calculating the initial rotor angle using the power flow equation.

```julia
    # Calculate currents from desired power
    Id0 = (P0 * Vd0 + Q0 * Vq0) / (Vd0^2 + Vq0^2)
    Iq0 = (P0 * Vq0 - Q0 * Vd0) / (Vd0^2 + Vq0^2)
```

These equations calculate the d-axis and q-axis currents required to produce the specified active and reactive power.

```julia
    # Calculate subtransient EMFs
    Eq_pp0 = Vd0 + machine.Ra * Id0 + machine.Xqpp * Iq0
    Ed_pp0 = Vq0 - machine.Ra * Iq0 - machine.Xdpp * Id0
    
    # Calculate transient EMFs
    Eq_p0 = Eq_pp0 + (machine.Xdp - machine.Xdpp) * Id0
    Ed_p0 = Ed_pp0 - (machine.Xqp - machine.Xqpp) * Iq0
```

These equations calculate the initial subtransient and transient EMFs based on the terminal conditions and machine parameters.

```julia
    # Calculate required field voltage
    Vf0 = Eq_p0 + (machine.Xd - machine.Xdp) * Id0
    
    # AVR equilibrium values
    Vr0 = Vf0 * avr.Ke
    Rf0 = avr.Kf * Vf0 / avr.Te
```

The field voltage (Vf0) is calculated to produce the required transient EMF. Then the regulator voltage (Vr0) and rate feedback (Rf0) are determined to ensure AVR equilibrium.

```julia
    # Governor values in equilibrium
    Pg0 = P0
    Pm_hp0 = turbine_gov.F_hp * P0
    Pm_lp0 = turbine_gov.F_lp * P0
```

The governor and turbine powers are set to match the desired active power output.

```julia
    return [δ0, 1.0, Eq_p0, Ed_p0, Eq_pp0, Ed_pp0, Vr0, Vf0, Rf0, Pg0, Pm_hp0, Pm_lp0]
end
```

The function returns the complete state vector with all 12 variables initialized.

## Running Simulations

```julia
function run_simulation(perturbations, tspan=(0.0, 10.0))
    # Set up parameters
    machine = default_machine_params()
    avr = default_avr_params()
    turbine_gov = default_turbine_gov_params()
    network = default_network_params()
    
    # Pack parameters for ODE solver
    p = (machine, avr, turbine_gov, network, perturbations)
    
    # Set up initial conditions
    u0 = init_steady_state(machine, avr, turbine_gov, network)
    
    # Verify steady state by checking derivatives
    du0 = zeros(12)
    synchronous_machine_dynamics!(du0, u0, p, 0.0)
    
    println("Initial derivatives:")
    for i in 1:length(du0)
        println("du[$i] = $(du0[i])")
    end
    
    if maximum(abs.(du0)) > 1e-3
        @warn "System not at equilibrium. Maximum derivative: $(maximum(abs.(du0)))"
    end
    
    # Define ODE problem
    prob = ODEProblem(synchronous_machine_dynamics!, u0, tspan, p)
    
    # Solve ODE problem
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
    
    return sol
end
```

This function:
1. Sets up parameters for all components
2. Initializes the state vector
3. Verifies the steady-state by checking initial derivatives
4. Creates and solves the ODE problem
5. Returns the solution

## Perturbation Studies

Three types of perturbations are implemented:

```julia
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
```

Each function:
1. Defines a specific perturbation (voltage reference change, load change, or both)
2. Runs the simulation with this perturbation
3. Analyzes and plots the results

## Analyzing Results

```julia
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
```

This extracts each state variable from the solution object for analysis.

```julia
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
        
        # Calculate currents using matrix solution
        A = [machine.Ra -machine.Xqpp; machine.Xdpp machine.Ra]
        b = [Eq_pp[i] - Vd; Ed_pp[i] + Vq]
        currents = A \ b
        Id = currents[1]
        Iq = currents[2]
        
        # Calculate electrical power
        Pe[i] = Vd*Id + Vq*Iq + machine.Ra*(Id^2 + Iq^2)
    end
```

This recalculates the terminal voltage and electrical power at each time step for plotting.

```julia
    # Plot results
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
    plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(800, 600), legend=true)
end
```

This creates six plots showing:
1. Rotor speed response
2. Rotor angle response
3. Terminal voltage response
4. Field voltage response
5. Electrical power response
6. Mechanical power components (high and low pressure)

## Main Function

```julia
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
```

The main function runs all three perturbation tests and saves the resulting plots to files.

## Troubleshooting

If the model shows instability:

1. **Check initialization**: Ensure initial derivatives are close to zero
2. **Increase damping**: Try setting a higher D value (1.0-4.0)
3. **Verify AVR parameters**: Ka might be too high causing voltage oscillations
4. **Add limiters**: Implement more aggressive limiting on derivatives
5. **Check network model**: Ensure the machine-network interface is correct

## Conclusion

This model provides a comprehensive representation of synchronous machine dynamics in power systems. It can be used to study machine responses to various disturbances, tune controllers, and analyze stability.

For further improvements, we might consider:
- Adding saturation effects
- Implementing more detailed excitation models
- Including power system stabilizers
- Extending to multi-machine systems