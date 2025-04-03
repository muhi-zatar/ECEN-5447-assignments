module Plotting

export plot_results, plot_voltages, plot_rotor, plot_machine_states, plot_line_currents

using Plots
using Plots.PlotMeasures
using LinearAlgebra
using ..ComponentIndex
using ..Transformations: ri_dq, dq_ri
using ..MachineComponents
using ..AVRComponents
using ..GovernorComponents
using ..NetworkComponents

"""
    plot_results(sol, indices::StateIndices, output_dir::String="results")

Generate all plots for simulation results and save them to the output directory.
"""
function plot_results(sol, indices::StateIndices, output_dir::String="results")
    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # Generate and save all plots
    plot_voltages(sol, indices, output_dir)
    plot_rotor(sol, indices, output_dir)
    plot_machine_states(sol, indices, output_dir)
    plot_line_currents(sol, indices, output_dir)
    plot_power(sol, indices, output_dir)
    plot_avr_states(sol, indices, output_dir)
    plot_governor_states(sol, indices, output_dir)
    
    return "Plots saved to $output_dir"
end

"""
    ri_dq_vector(d_values, q_values)

Convert vectors of d and q values to magnitude and angle.
Always uses angle 0.0 for network reference frame.
"""
function ri_dq_vector(d_values, q_values)
    # Map the network D and Q values onto real and imaginary axes
    # Always use angle 0.0 for network reference frame
    RI_values = map((v_d, v_q) -> dq_ri(0.0) * [v_d; v_q], d_values, q_values)

    # Take the magnitude
    mag = map(V -> hypot(V[1], V[2]), RI_values)

    return RI_values, mag
end

"""
    compute_S_vector(v_RI_values, i_RI_values)

Compute apparent power from voltage and current vectors.
"""
function compute_S_vector(v_RI_values, i_RI_values)
    # Compute apparent power based on vectors of real and imaginary voltage and current
    S_values = map((V, I) -> [
            V[1] * I[1] + V[2] * I[2];  # P = V_R * I_R + V_I * I_I
            V[2] * I[1] - V[1] * I[2]   # Q = V_I * I_R - V_R * I_I
        ], v_RI_values, i_RI_values)

    return S_values
end

"""
    plot_voltages(sol, indices::StateIndices, output_dir::String)

Generate bus voltage magnitude plots.
"""
function plot_voltages(sol, indices::StateIndices, output_dir::String)
    t = sol.t
    
    # Extract voltage components
    v_1_d_values = [sol[indices.network[:V_1_D_IDX], i] for i in 1:length(t)]
    v_1_q_values = [sol[indices.network[:V_1_Q_IDX], i] for i in 1:length(t)]
    v_2_d_values = [sol[indices.network[:V_2_D_IDX], i] for i in 1:length(t)]
    v_2_q_values = [sol[indices.network[:V_2_Q_IDX], i] for i in 1:length(t)]
    v_3_d_values = [sol[indices.network[:V_3_D_IDX], i] for i in 1:length(t)]
    v_3_q_values = [sol[indices.network[:V_3_Q_IDX], i] for i in 1:length(t)]
    
    # Transform to RI and calculate magnitudes
    _, v_1_magnitude = ri_dq_vector(v_1_d_values, v_1_q_values)
    _, v_2_magnitude = ri_dq_vector(v_2_d_values, v_2_q_values)
    _, v_3_magnitude = ri_dq_vector(v_3_d_values, v_3_q_values)
    
    # Create plot
    p = plot(t, v_1_magnitude, 
        label="Infinite Bus", 
        linewidth=2, 
        title="Bus Voltage Magnitudes", 
        xlabel="Time (s)",
        ylabel="Voltage (pu)",
        left_margin=10mm)
    plot!(p, t, v_2_magnitude, label="Machine Bus", linewidth=2)
    plot!(p, t, v_3_magnitude, label="Load Bus", linewidth=2)
    
    # Save plot
    savefig(p, joinpath(output_dir, "bus_voltages.png"))
    
    return p
end

"""
    plot_rotor(sol, indices::StateIndices, output_dir::String)

Generate rotor angle and speed plots.
"""
function plot_rotor(sol, indices::StateIndices, output_dir::String)
    t = sol.t
    
    # Extract rotor states
    delta_values = [sol[indices.machine[:DELTA], i] for i in 1:length(t)]
    omega_values = [sol[indices.machine[:OMEGA], i] for i in 1:length(t)]
    
    # Create angle plot
    p1 = plot(t, delta_values,
        label="Rotor angle (δ)", 
        title="Rotor Angle", 
        linewidth=2, 
        xlabel="Time (s)",
        ylabel="Angle (rad)",
        left_margin=10mm)
    
    # Create speed plot
    p2 = plot(t, omega_values,
        label="Rotor speed (ω)", 
        title="Rotor Speed", 
        linewidth=2, 
        xlabel="Time (s)",
        ylabel="Speed (pu)",
        left_margin=10mm)
    
    # Save plots
    savefig(p1, joinpath(output_dir, "rotor_angle.png"))
    savefig(p2, joinpath(output_dir, "rotor_speed.png"))
    
    # Combined plot
    p_combined = plot(p1, p2, layout=(2,1), size=(800, 600))
    savefig(p_combined, joinpath(output_dir, "rotor_states.png"))
    
    return p_combined
end

"""
    plot_machine_states(sol, indices::StateIndices, output_dir::String)

Generate machine flux and other state plots.
"""
function plot_machine_states(sol, indices::StateIndices, output_dir::String)
    t = sol.t
    
    # Extract flux states
    eq_p_values = [sol[indices.machine[:EQ_P], i] for i in 1:length(t)]
    ed_p_values = [sol[indices.machine[:ED_P], i] for i in 1:length(t)]
    psi_d_pp_values = [sol[indices.machine[:PSI_D_PP], i] for i in 1:length(t)]
    psi_q_pp_values = [sol[indices.machine[:PSI_Q_PP], i] for i in 1:length(t)]
    psi_d_values = [sol[indices.machine[:PSI_D], i] for i in 1:length(t)]
    psi_q_values = [sol[indices.machine[:PSI_Q], i] for i in 1:length(t)]
    
    # Create flux plots
    p = plot(t, eq_p_values,
        label="eq'", 
        title="Machine Fluxes", 
        linewidth=2, 
        xlabel="Time (s)",
        ylabel="Flux (pu)",
        left_margin=10mm)
    plot!(p, t, ed_p_values, label="ed'", linewidth=2)
    plot!(p, t, psi_d_pp_values, label="ψd''", linewidth=2)
    plot!(p, t, psi_q_pp_values, label="ψq''", linewidth=2)
    plot!(p, t, psi_d_values, label="ψd", linewidth=2)
    plot!(p, t, psi_q_values, label="ψq", linewidth=2)
    
    # Save plot
    savefig(p, joinpath(output_dir, "machine_fluxes.png"))
    
    return p
end

"""
    plot_line_currents(sol, indices::StateIndices, output_dir::String)

Generate line current magnitude plots.
"""
function plot_line_currents(sol, indices::StateIndices, output_dir::String)
    t = sol.t
    
    # Extract line current components
    i_12_d_values = [sol[indices.network[:I_12_D_IDX], i] for i in 1:length(t)]
    i_12_q_values = [sol[indices.network[:I_12_Q_IDX], i] for i in 1:length(t)]
    i_23_d_values = [sol[indices.network[:I_23_D_IDX], i] for i in 1:length(t)]
    i_23_q_values = [sol[indices.network[:I_23_Q_IDX], i] for i in 1:length(t)]
    i_13_d_values = [sol[indices.network[:I_13_D_IDX], i] for i in 1:length(t)]
    i_13_q_values = [sol[indices.network[:I_13_Q_IDX], i] for i in 1:length(t)]
    
    # Transform to RI and calculate magnitudes
    _, i_12_mag = ri_dq_vector(i_12_d_values, i_12_q_values)
    _, i_23_mag = ri_dq_vector(i_23_d_values, i_23_q_values)
    _, i_13_mag = ri_dq_vector(i_13_d_values, i_13_q_values)
    
    # Create plot
    p = plot(t, i_12_mag, 
        label="Line 1-2", 
        title="Line Currents", 
        linewidth=2, 
        xlabel="Time (s)",
        ylabel="Current (pu)",
        left_margin=10mm)
    plot!(p, t, i_23_mag, label="Line 2-3", linewidth=2)
    plot!(p, t, i_13_mag, label="Line 1-3", linewidth=2)
    
    # Save plot
    savefig(p, joinpath(output_dir, "line_currents.png"))
    
    return p
end

"""
    plot_power(sol, indices::StateIndices, output_dir::String)

Generate bus power plots.
"""
function plot_power(sol, indices::StateIndices, output_dir::String)
    t = sol.t
    
    # Extract voltage and current components for buses 1 and 3
    v_1_d_values = [sol[indices.network[:V_1_D_IDX], i] for i in 1:length(t)]
    v_1_q_values = [sol[indices.network[:V_1_Q_IDX], i] for i in 1:length(t)]
    v_3_d_values = [sol[indices.network[:V_3_D_IDX], i] for i in 1:length(t)]
    v_3_q_values = [sol[indices.network[:V_3_Q_IDX], i] for i in 1:length(t)]
    
    i_1_d_values = [sol[indices.network[:I_1_D_IDX], i] for i in 1:length(t)]
    i_1_q_values = [sol[indices.network[:I_1_Q_IDX], i] for i in 1:length(t)]
    i_3_d_values = [sol[indices.network[:I_3_D_IDX], i] for i in 1:length(t)]
    i_3_q_values = [sol[indices.network[:I_3_Q_IDX], i] for i in 1:length(t)]
    
    # Transform to RI
    v_1_RI_values, _ = ri_dq_vector(v_1_d_values, v_1_q_values)
    v_3_RI_values, _ = ri_dq_vector(v_3_d_values, v_3_q_values)
    i_1_RI_values, _ = ri_dq_vector(i_1_d_values, i_1_q_values)
    i_3_RI_values, _ = ri_dq_vector(i_3_d_values, i_3_q_values)
    
    # Compute power
    S_1_values = compute_S_vector(v_1_RI_values, i_1_RI_values)
    S_3_values = compute_S_vector(v_3_RI_values, i_3_RI_values)
    
    # Create plot
    p = plot(t, first.(S_1_values), 
        label="Inf. Bus (P)", 
        title="Bus Power", 
        linewidth=2, 
        xlabel="Time (s)",
        ylabel="Power (pu)",
        left_margin=10mm)
    plot!(p, t, last.(S_1_values), label="Inf. Bus (Q)", linewidth=2)
    plot!(p, t, first.(S_3_values), label="Load (P)", linewidth=2)
    plot!(p, t, last.(S_3_values), label="Load (Q)", linewidth=2)
    
    # Save plot
    savefig(p, joinpath(output_dir, "powers.png"))
    
    return p
end

"""
    plot_avr_states(sol, indices::StateIndices, output_dir::String)

Generate AVR state plots.
"""
function plot_avr_states(sol, indices::StateIndices, output_dir::String)
    t = sol.t
    
    # Extract AVR states
    efd_values = [sol[indices.avr[:EFD_IDX], i] for i in 1:length(t)]
    vs_values = [sol[indices.avr[:VS_IDX], i] for i in 1:length(t)]
    vll_values = [sol[indices.avr[:VLL_IDX], i] for i in 1:length(t)]
    vf_values = [sol[indices.avr[:VF_IDX], i] for i in 1:length(t)]
    
    # Create plot
    p = plot(t, efd_values, 
        label="Field Voltage", 
        title="AVR States", 
        linewidth=2, 
        xlabel="Time (s)",
        ylabel="Value (pu)",
        left_margin=10mm)
    plot!(p, t, vs_values, label="Sensed Terminal Voltage", linewidth=2)
    plot!(p, t, vll_values, label="Lead-Lag State", linewidth=2)
    plot!(p, t, vf_values, label="Feedback State", linewidth=2)
    
    # Save plot
    savefig(p, joinpath(output_dir, "avr_states.png"))
    
    return p
end

"""
    plot_governor_states(sol, indices::StateIndices, output_dir::String)

Generate governor state plots.
"""
function plot_governor_states(sol, indices::StateIndices, output_dir::String)
    t = sol.t
    
    # Extract governor states
    fv_values = [sol[indices.governor[:FV_IDX], i] for i in 1:length(t)]
    ff_values = [sol[indices.governor[:FF_IDX], i] for i in 1:length(t)]
    et_values = [sol[indices.governor[:ET_IDX], i] for i in 1:length(t)]
    
    # Create plot
    p = plot(t, fv_values, 
        label="Fuel Value", 
        title="Governor States", 
        linewidth=2, 
        xlabel="Time (s)",
        ylabel="Value (pu)",
        left_margin=10mm)
    plot!(p, t, ff_values, label="Fuel Flow", linewidth=2)
    plot!(p, t, et_values, label="Exhaust Temp", linewidth=2)
    
    # Save plot
    savefig(p, joinpath(output_dir, "gov_states.png"))
    
    return p
end

end # module