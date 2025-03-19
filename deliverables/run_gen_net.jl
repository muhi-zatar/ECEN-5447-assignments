cd(@__DIR__)
using Pkg
Pkg.activate("../.")

using DifferentialEquations
using Plots
using LinearAlgebra
using Sundials

include("NetworkModel.jl")
include("SauerPaiMachineModel.jl")
include("PowerFlowWrapper.jl")

using .NetworkModel
using .SauerPaiMachineModel
using .PowerFlowWrapper


const DELTA = 1
const OMEGA = 2
const EQ_P = 3
const ED_P = 4
const PSI_D_PP = 5
const PSI_Q_PP = 6

mutable struct MachineNetworkParams
    network::ThreeBusNetwork
    machine::SauerPaiMachine
    Vf::Float64                 
    τm::Float64                  
    network_idx::UnitRange{Int64}
    machine_idx::UnitRange{Int64}
    θ::Float64
    i_2_r::Float64
    i_2_i::Float64
end

function run_machine_network(network_file)
    V_sol, θ_sol, P_sol, Q_sol = get_powerflow_results(network_file)

    V_mag = V_sol[BUS_MACHINE_MODEL]
    V_angle = θ_sol[BUS_MACHINE_MODEL]
    P = P_sol[BUS_MACHINE_MODEL]
    Q = Q_sol[BUS_MACHINE_MODEL]

    V_terminal = V_mag * exp(im * V_angle)
    println("Initial power flow results:")
    println("V_terminal = $V_mag p.u ∠ $V_angle rad")
    println("P = $P pu, Q = $Q pu")

    network = ThreeBusNetwork()
    
    machine = SauerPaiMachine(
        base_power=100.0,
        system_base_power=100.0,
        system_base_frequency=60.0
    )

    network_states, i_2_r_init, i_2_i_init = initialize_network(network, V_sol, θ_sol, P_sol, Q_sol)
    machine_states, Vf_init, τ_m_init = initialize_machine(machine, V_terminal, V_angle, P, Q)
    
    println("\nInitial States:")
    println("Machine states:")
    println("  Delta (rotor angle): $(machine_states[DELTA])")
    println("  Omega (rotor speed): $(machine_states[OMEGA])")
    println("  EQ_P: $(machine_states[EQ_P])")
    println("  ED_P: $(machine_states[ED_P])")
    println("  PSI_D_PP: $(machine_states[PSI_D_PP])")
    println("  PSI_Q_PP: $(machine_states[PSI_Q_PP])")
    
    println("Initial field voltage (Vf): $Vf_init")
    println("Initial mechanical torque (τm): $τ_m_init")
    println("Initial current injection: real=$i_2_r_init, imag=$i_2_i_init")

    states = vcat(network_states, machine_states)
    
    # Define state 
    network_idx = 1:22
    machine_idx = 23:28

    p = MachineNetworkParams(
        network,
        machine,
        Vf_init,     
        τ_m_init,    
        network_idx,
        machine_idx,
        V_angle,      
        i_2_r_init,   
        i_2_i_init    
    )

    M_system = zeros(Float64, length(states))
    
    if isa(network.M, Vector)
        M_system[network_idx] .= network.M
    else
        for i in 1:length(network_idx)
            M_system[network_idx[i]] = network.M[i, i]
        end
    end
    
    M_system[machine_idx] .= 1.0
    
    mass_matrix = Diagonal(M_system)

    function machine_network_dynamics!(du, u, params, t, mass_matrix=M_system)
        network_states = @view u[params.network_idx]
        machine_states = @view u[params.machine_idx]
        
        du_network = similar(network_states, length(params.network_idx))
        du_machine = similar(machine_states, length(params.machine_idx))
        
        network_states_f64 = convert.(Float64, network_states)
        
        V_terminal = update_network_states!(
            network_states_f64,
            du_network,
            params.i_2_r,
            params.i_2_i,
            params.θ,
            params.network
        )
        
        V_mag, I_RI, θ, ω = update_machine_states!(
            machine_states,
            du_machine,
            V_terminal,
            params.Vf,
            params.τm,
            params.machine
        )
        
        i_2_r = real(I_RI)
        i_2_i = imag(I_RI)
        
        for i in 1:length(params.network_idx)
            du[params.network_idx[i]] = du_network[i]
        end
        du[params.machine_idx] .= du_machine
        
        params.θ = θ
        params.i_2_r = i_2_r
        params.i_2_i = i_2_i
        
        if abs(t - round(t)) < 0.001
            println("t=$t: δ=$(machine_states[DELTA]), ω=$ω, τm=$(params.τm), Vf=$(params.Vf), V_mag=$V_mag")
        end
    end


    tspan = (0.0, 10.0)  
    prob = ODEProblem(machine_network_dynamics!, states, tspan, p)
    
    sol = solve(prob, Rodas5(autodiff=false), reltol=1e-6, abstol=1e-6)

    t = sol.t
    
    p1 = plot(t, [sol[machine_idx[1], i] for i in 1:length(t)], 
              label="Rotor angle (δ)", title="Machine States", linewidth=2)
    plot!(p1, t, [sol[machine_idx[2], i] for i in 1:length(t)], 
          label="Rotor speed (ω)", linewidth=2)
    
    p2 = plot(t, [sol[machine_idx[3], i] for i in 1:length(t)], 
              label="eq'", title="Machine Fluxes", linewidth=2)
    plot!(p2, t, [sol[machine_idx[4], i] for i in 1:length(t)], 
          label="ed'", linewidth=2)
    plot!(p2, t, [sol[machine_idx[5], i] for i in 1:length(t)], 
          label="ψd''", linewidth=2)
    plot!(p2, t, [sol[machine_idx[6], i] for i in 1:length(t)], 
          label="ψq''", linewidth=2)
    
    p3 = plot(t, [sol[network_idx[9], i] for i in 1:length(t)], 
              title="Terminal Voltage (Machine Bus)", label="Vd", linewidth=2)
    plot!(p3, t, [sol[network_idx[10], i] for i in 1:length(t)], 
          label="Vq", linewidth=2)
    
    p4 = plot(t, [p.i_2_r for i in 1:length(t)], 
              title="Current Injection (Machine Bus)", label="Id", linewidth=2)
    plot!(p4, t, [p.i_2_i for i in 1:length(t)], 
          label="Iq", linewidth=2)
    
    p_combined = plot(p1, p2, p3, p4, layout=(4,1), size=(800, 1600))
    savefig(p_combined, "machine_network_results.png")
    
    return sol
end

sol = run_machine_network("../data/ThreeBusMultiLoad.raw")  # Adjust path as needed

println("Simulation complete!")