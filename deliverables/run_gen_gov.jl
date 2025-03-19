cd(@__DIR__)
using Pkg
Pkg.activate("../.")

using DifferentialEquations
using Plots
using LinearAlgebra

include("SauerPaiMachineModel.jl")
include("GasTGModel.jl")

using .SauerPaiMachineModel
using .GasTGModel


const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6

const FV_IDX = 1                # Fuel value
const FF_IDX = 2                # Fuel flow
const ET_IDX = 3                # Exhaust temp

struct MachineGovParams
    machine::SauerPaiMachine
    governor::GasTG
    V_terminal::Complex{Float64}  # Fixed terminal voltage
    Vf::Float64                   # Fixed field voltage
    machine_idx::UnitRange{Int64}
    gov_idx::UnitRange{Int64}
end

function run_machine_governor()
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

    governor = GasTG(
        P_ref=P 
    )

    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)
    ω_init = 1.0
    governor_states = initialize_gov(governor, τ_m_init, ω_init)
    
    println("\nInitial States:")
    println("Machine states:")
    println("  Delta (rotor angle): $(machine_states[DELTA])")
    println("  Omega (rotor speed): $(machine_states[OMEGA])")
    println("  EQ_P: $(machine_states[EQ_P])")
    println("  ED_P: $(machine_states[ED_P])")
    println("  PSI_D_PP: $(machine_states[PSI_D_PP])")
    println("  PSI_Q_PP: $(machine_states[PSI_Q_PP])")
    
    println("Governor states:")
    println("  Fuel Value: $(governor_states[FV_IDX])")
    println("  Fuel Flow: $(governor_states[FF_IDX])")
    println("  Exhaust Temp: $(governor_states[ET_IDX])")
    
    println("Initial field voltage (Vf): $Vf_init")
    println("Initial mechanical torque (τm): $τ_m_init")
    states = vcat(machine_states, governor_states)
    
    machine_idx = 1:6
    gov_idx = 7:9

    p = MachineGovParams(
        machine,
        governor,
        V_terminal,
        Vf_init,      
        machine_idx,
        gov_idx
    )

    function machine_gov_dynamics!(du, u, params, t)
        machine_states = @view u[params.machine_idx]
        gov_states = @view u[params.gov_idx]
        
        du_machine = zeros(Float64, length(params.machine_idx))
        du_gov = zeros(Float64, length(params.gov_idx))
        
        ω = machine_states[OMEGA]
        
        τm = update_gov_states!(
            gov_states,
            du_gov,
            ω,
            params.governor
        )
        
        V_mag, I_RI, θ, ω_updated = update_machine_states!(
            machine_states,
            du_machine,
            params.V_terminal,
            params.Vf,
            τm,
            params.machine
        )
        
        du[params.machine_idx] .= du_machine
        du[params.gov_idx] .= du_gov
        
        if abs(t - round(t)) < 0.001
            println("t=$t: δ=$(machine_states[DELTA]), ω=$ω, τm=$τm, Vf=$(params.Vf), V_mag=$V_mag")
        end
    end

    tspan = (0.0, 10.0)
    prob = ODEProblem(machine_gov_dynamics!, states, tspan, p)
    
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

    t = sol.t
    
    delta_values = [sol[machine_idx[DELTA], i] for i in 1:length(t)]
    omega_values = [sol[machine_idx[OMEGA], i] for i in 1:length(t)]
    
    τm_values = []
    τe_values = []
    
    for i in 1:length(t)
        gov_states_i = [sol[gov_idx[j], i] for j in 1:length(gov_idx)]
        
        
        τm = sol[gov_idx[FV_IDX], i]  # Fuel value is approximately mechanical torque
        push!(τm_values, τm)
        
        # Simplified approximation of electrical torque
        delta = sol[machine_idx[DELTA], i]
        eq_p = sol[machine_idx[EQ_P], i]
        ed_p = sol[machine_idx[ED_P], i]
        
        # This is a simplified approximation
        τe = eq_p * sin(delta) - ed_p * cos(delta)
        push!(τe_values, τe)
    end
    
    p1 = plot(t, delta_values, 
              label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, omega_values, 
          label="Rotor speed (ω)", linewidth=2)
    
    p2 = plot(t, [sol[machine_idx[EQ_P], i] for i in 1:length(t)], 
              label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[machine_idx[ED_P], i] for i in 1:length(t)], 
          label="ed'", linewidth=2)
    plot!(p2, t, [sol[machine_idx[PSI_D_PP], i] for i in 1:length(t)], 
          label="ψd''", linewidth=2)
    plot!(p2, t, [sol[machine_idx[PSI_Q_PP], i] for i in 1:length(t)], 
          label="ψq''", linewidth=2)
    
    p3 = plot(t, [sol[gov_idx[FV_IDX], i] for i in 1:length(t)], 
              label="Fuel Value", title="Governor States", linewidth=2)
    plot!(p3, t, [sol[gov_idx[FF_IDX], i] for i in 1:length(t)], 
          label="Fuel Flow", linewidth=2)
    plot!(p3, t, [sol[gov_idx[ET_IDX], i] for i in 1:length(t)], 
          label="Exhaust Temp", linewidth=2)
    
    p4 = plot(t, τe_values, 
              label="Electrical Torque (approx)", title="Machine Torques", linewidth=2)
    plot!(p4, t, τm_values, 
          label="Mechanical Torque", linewidth=2)
    
    p_combined = plot(p1, p2, p3, p4, layout=(4,1), size=(800, 1600))
    savefig(p_combined, "machine_governor_results.png")
    
    return sol
end

sol = run_machine_governor()

println("Simulation complete!")