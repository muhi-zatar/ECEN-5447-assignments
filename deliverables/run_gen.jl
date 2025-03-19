cd(@__DIR__)
using Pkg
Pkg.activate("../.")

using DifferentialEquations
using Plots
using LinearAlgebra

include("SauerPaiMachineModel.jl")

using .SauerPaiMachineModel

const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6

struct MachineParams
    machine::SauerPaiMachine
    V_terminal::Complex{Float64}  # Fixed terminal voltage
    Vf::Float64                  # Fixed field voltage
    τm::Float64                  # Fixed mechanical torque
end

function run_machine_only()
    V_mag = 1.0142        
    V_angle = -0.007843738409384566      
    P = 0.8          
    Q = 0.10751144041475173
    
    V_terminal = V_mag * exp(im * V_angle)
    println("Initial conditions:")
    println("V_terminal = $V_mag p.u ∠ $V_angle rad")
    println("P = $P pu, Q = $Q pu")

    machine = SauerPaiMachine(
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)
    
    println("\nInitial Machine States:")
    println("Delta (rotor angle): $(machine_states[DELTA])")
    println("Omega (rotor speed): $(machine_states[OMEGA])")
    println("EQ_P: $(machine_states[EQ_P])")
    println("ED_P: $(machine_states[ED_P])")
    println("PSI_D_PP: $(machine_states[PSI_D_PP])")
    println("PSI_Q_PP: $(machine_states[PSI_Q_PP])")
    println("Initial field voltage (Vf): $Vf_init")
    println("Initial mechanical torque (τm): $τ_m_init")

    # Create params for ODE solver
    p = MachineParams(
        machine,
        V_terminal, 
        Vf_init,
        τ_m_init
    )

    function machine_dynamics!(du, u, params, t)
        machine_states = u
        
        du_machine = zeros(Float64, 6)
        
        V_mag, I_RI, θ, ω = update_machine_states!(
            machine_states,
            du_machine,
            params.V_terminal,
            params.Vf,
            params.τm,
            params.machine
        )
        
        du .= du_machine
        
        if abs(t - round(t)) < 0.001
            println("t=$t: δ=$(machine_states[DELTA]), ω=$ω, τm=$(params.τm), Vf=$(params.Vf), V_mag=$V_mag")
        end
    end

    tspan = (0.0, 20.0) 
    prob = ODEProblem(machine_dynamics!, machine_states, tspan, p)
    
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

    t = sol.t
    
    p1 = plot(t, [sol[DELTA, i] for i in 1:length(t)], 
              label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, [sol[OMEGA, i] for i in 1:length(t)], 
          label="Rotor speed (ω)", linewidth=2)
    
    p2 = plot(t, [sol[EQ_P, i] for i in 1:length(t)], 
              label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[ED_P, i] for i in 1:length(t)], 
          label="ed'", linewidth=2)
    plot!(p2, t, [sol[PSI_D_PP, i] for i in 1:length(t)], 
          label="ψd''", linewidth=2)
    plot!(p2, t, [sol[PSI_Q_PP, i] for i in 1:length(t)], 
          label="ψq''", linewidth=2)
    
    τe_values = []
    for i in 1:length(t)
        # This calculation requires reconstructing the electrical torque from states
        delta = sol[DELTA, i]
        eq_p = sol[EQ_P, i]
        ed_p = sol[ED_P, i]
        
        # This is a approximation
        τe = eq_p * sin(delta) - ed_p * cos(delta)
        push!(τe_values, τe)
    end
    
    p3 = plot(t, τe_values, 
              label="Electrical Torque (approx)", title="Machine Torques", linewidth=2)
    plot!(p3, t, fill(p.τm, length(t)), 
          label="Mechanical Torque (constant)", linewidth=2)
    
    p_combined = plot(p1, p2, p3, layout=(3,1), size=(800, 1200))
    savefig(p_combined, "machine_only_results.png")
    
    return sol
end

sol = run_machine_only()

println("Simulation complete!")