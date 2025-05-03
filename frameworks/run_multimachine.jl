cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.resolve() # make sure Manifest matches Project
Pkg.instantiate() # install missing dependencies

# Importing necessary Modules
using PowerSystemCaseBuilder
using PowerSimulationsDynamics
using PowerNetworkMatrices
using Plots
using Plots.PlotMeasures
using Sundials
using InfrastructureSystems
using PowerSystems
using Logging
const PSY = PowerSystems

# Build the system frp, PSID systems, this should be similar to ours.
sys = build_system(PSIDTestSystems, "psid_test_threebus_multimachine_dynlines")

# Define a helper function for plotting
function plot_stuff(results, title=nothing)
    # Bus voltages
    V_mag_series_SM = get_voltage_magnitude_series(results, 101)
    V_mag_series_CV = get_voltage_magnitude_series(results, 102)
    V_mag_series_L = get_voltage_magnitude_series(results, 103)

    p1 = plot(V_mag_series_SM, label="SG Bus", title="Bus Voltages", ylim=(0.85, 1.15), linewidth=2)
    plot!(p1, V_mag_series_CV, label="CV Bus", linewidth=2)
    plot!(p1, V_mag_series_L, label="Load Bus", linewidth=2)
    ylabel!(p1, "Voltage [p.u.]")
    xlabel!(p1, "Time [s]")
    if title !== nothing
        savefig(p1, "Voltage_$title.png")
    else
        savefig(p1, "Voltage.png")
    end

    # Rotor angle
    rotor_angle_series = get_state_series(results, ("generator-101-1", :δ))
    pll_angle_series = get_state_series(results, ("generator-102-1", :θ_pll))
    oc_angle_series = get_state_series(results, ("generator-102-1", :θ_oc))

    p2 = plot(rotor_angle_series[1], rad2deg.(rotor_angle_series[2]), label="SG Bus", title="Rotor Angles", ylim=(5, 20), linewidth=2)
    plot!(p2, oc_angle_series[1], rad2deg.(oc_angle_series[2]), label="OLC", linewidth=2)
    plot!(p2, pll_angle_series[1], rad2deg.(pll_angle_series[2]), label="PLL", linewidth=2)

    ylabel!(p2, "δ [degree]")
    xlabel!(p2, "Time [s]")
    if title !== nothing
        savefig(p2, "rotor_angle_$title.png")
    else
        savefig(p2, "rotor_angle.png")
    end

    # Frequency
    rotor_speed_series = get_state_series(results, ("generator-101-1", :ω))
    outerloop_power_series = get_state_series(results, ("generator-102-1", :p_oc))
    outerloop_ω_ref = 1.0
    outerloop_Rp = 0.05
    outerloop_P_ref = results.setpoints["generator-102-1"]["P_ref"]
    ω_olc_series = (outerloop_ω_ref .+ outerloop_Rp) .* (outerloop_P_ref .- outerloop_power_series[2])

    p3 = plot(rotor_speed_series[1], (rotor_speed_series[2] .* 60), label="SG Bus", title="Rotor Speed", ylim=(59.9, 60.1), linewidth=2)
    plot!(p3, outerloop_power_series[1], (ω_olc_series .* 60), label="OLC", linewidth=2)
    ylabel!(p3, "ω [Hz]")
    xlabel!(p3, "Time [s]")
    if title !== nothing
        savefig(p3, "rotor_speed_$title.png")
    else
        savefig(p3, "rotor_speed.png")
    end

    # Active power
    active_power_series_SM = get_activepower_series(results, "generator-101-1")
    active_power_series_CV = get_activepower_series(results, "generator-102-1")
    active_power_series_L = get_activepower_series(results, "load1031")

    p4 = plot(active_power_series_SM[1], ((active_power_series_SM[2] .- active_power_series_SM[2][1]) .* 100), label="SG Bus", title="Bus Active Power Generation/Consumption", ylim=(-50, 50), linewidth=2)
    plot!(p4, active_power_series_CV[1], ((active_power_series_CV[2] .- active_power_series_CV[2][1]) .* 100), label="CV Bus", linewidth=2)
    plot!(p4, active_power_series_L[1], ((active_power_series_L[2] .- active_power_series_L[2][1]) .* 100), label="Load Bus", linewidth=2)
    ylabel!(p4, "Active Power Deviation [MW]")
    xlabel!(p4, "Time [s]")
    if title !== nothing
        savefig(p4, "active_power_$title.png")
    else
        savefig(p4, "active_power.png")
    end
    return
end

##### STEP 1: MODIFY THE NETWORK TO MATCH OUR MODEL #####

# See which buses have generators
machine_gen = get_component(ThermalStandard, sys, "generator-101-1")        # Synchronous Machine goes on Bus 1
converter_gen = get_component(ThermalStandard, sys, "generator-102-1")      # Inverter goes on Bus 2

# Remove dynamic injectors (to prepare for replacement by our model)
remove_component!(sys, get_dynamic_injector(converter_gen))
remove_component!(sys, get_dynamic_injector(machine_gen))

# Define the slack bus
slack_bus = first(get_components(x -> get_bustype(x) == ACBusTypes.REF, Bus, sys))

# Ensure that machine goes on slack bus
if machine_gen.bus != slack_bus
    @warn("Placing the machine on a different bus than the slack bus will cause the ResidualModel build to fail!")
end

### Add a dynamic generator

# Define SauerPaiMachine parameters
sp_machine = SauerPaiMachine(
    0.002, # R
    1.79, # Xd
    1.71, # Xq
    0.169, # Xd_p
    0.228, # Xq_p
    0.135, # Xd_pp
    0.2, # Xq_pp
    0.13, # Xl
    4.3, # Td0_p
    0.85, # Tq0_p 
    0.032, # Td0_pp
    0.05, # Tq0_pp
)

# Define Shaft with single mass
shaft_model = SingleMass(
    H=3.148,
    D=2.0,
)

# Define GasTG governer
governor = GasTG(
    0.05,    # R
    0.2,     # T1
    0.2,     # T2
    2.0,    # T3
    1.0,     # AT
    2.5,     # Kt
    (0.01, 1.1), # Vlim
    0.0    # D_turb
)

# Define EXST1 AVR
avr = PSY.EXST1(
    0.01, # Tr
    (-5.0, 5.0), # Vi_lim
    10.0, # Tc
    20.0, # Tb
    200.0, # Ka
    0.1, # Ta
    (0.0, 6.0), # Vr_limit
    0.0, # Kc
    0.0, # Kf
    0.1, # Tf
    1.0, # V_ref
)

# Define dynamic generator wrapper
dg1 = DynamicGenerator(
    name=get_name(machine_gen),
    ω_ref=1.0,
    machine=sp_machine,
    shaft=shaft_model,
    avr=avr,
    prime_mover=governor,
    pss=PSY.PSSFixed(V_pss=0.0),
    base_power=get_base_power(machine_gen),
    # InfrastructureSystems.InfrastructureSystemsInternal()
)

# Add our generator to the system
add_component!(sys, dg1, machine_gen)

### Add a dynamic inverter

# Define converter as an AverageConverter
converter_high_power() = AverageConverter(
    rated_voltage=138.0,
    rated_current=100.0
)

# Define Outer Control as composition of Real + Reactive Power droop
outer_control() = OuterControl(
    ActivePowerDroop(Rp=0.05, ωz=2π * 60),
    ReactivePowerDroop(kq=2.0, ωf=0.1322 * π * 50.0)     # TODO: Understand Kpq/Kiq distinction
)

# Define Inner Control as a Voltage+Current Controller with Virtual Impedance
inner_control() = VoltageModeControl(
    kpv=0.59,
    kiv=736.0,
    kffv=0.0,
    rv=0.0,
    lv=0.2,
    kpc=1.27,
    kic=14.3,
    kffi=0.0,
    ωad=50.0,
    kad=0.2
)

# Define DC Source as a FixedSource:
dc_source_lv() = FixedDCSource(voltage=600.0)

# Define a Frequency Estimator as a PLL
pll() = ReducedOrderPLL(
    ω_lp=1.322 * π * 50,
    kp_pll=2.0,
    ki_pll=20.0,
)

# Define an LCL filter:
filt() = LCLFilter(lf=0.08, rf=0.003, cf=0.074, lg=0.2, rg=0.01)

# Define dynamic inverter wrapper
case_inv = DynamicInverter(
    get_name(converter_gen),
    1.0,
    converter_high_power(),
    outer_control(),
    inner_control(),
    dc_source_lv(),
    pll(),
    filt(),
)

# Add our dynamic inverter to the system
add_component!(sys, case_inv, converter_gen)

##### STEP 2: RUN A DYNAMIC SIMULATION IN STEADY-STATE #####
# Define Simulation
tspan = (0.0, 30.0)
sim = Simulation(
    ResidualModel,
    sys,
    mktempdir(),
    tspan,
    console_level=Logging.Info
)

# Print the initial values
show_states_initial_value(sim)

# execute the simulation
execute!(
    sim, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax=0.02, #Maximum step size
)

results = read_results(sim)

### Plot some things
plot_stuff(results, "steady-state")

##### STEP 3: RUN A SIMULATION FOR A LINE TRIP #####

#Make a copy of the original system
new_sys = deepcopy(sys);

#Remove Line "BUS 1-BUS 2-i_1"
remove_component!(DynamicBranch, new_sys, "BUS 1-BUS 2-i_1")

# Obtain the new Ybus
Ybus_fault_dyn = Ybus(new_sys).data

Ybus_change_dyn = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault_dyn, #New YBus
)

# Define Simulation
sim_dyn = Simulation(
    ResidualModel, #Type of model used
    sys, #system
    pwd(), #folder to output results
    (0.0, 5.0), #time span
    Ybus_change_dyn, #Type of perturbation
    ; console_level=Logging.Debug
)

#Run the simulation
execute!(
    sim_dyn, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax=0.02, #Maximum step size
)

line_trip_results = read_results(sim_dyn)

# Plot some things
plot_stuff(line_trip_results, "line_trip")

##### STEP 4: RUN A SIMULATION OF A LOAD INCREASE #####
# Get the device
l_device = get_component(ElectricLoad, sys, "load1031")

# Increase the load by decreasing the reference value
load_increase = LoadChange(1.0, l_device, :P_ref, (1.8 * 0.85))

# Define Simulation
sim_dyn_load_increase = Simulation(
    ResidualModel, #Type of model used
    sys, #system
    pwd(), #folder to output results
    (0.0, 30.0), #time span
    load_increase, #Type of perturbation
)

#Run the simulation
execute!(
    sim_dyn_load_increase, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax=0.02, #Maximum step size
)

load_increase_results = read_results(sim_dyn_load_increase)

# Plot results
plot_stuff(load_increase_results, "load_increase")


##### STEP 5: RUN A SIMULATION OF A LOAD DECREASE #####
# Get the device
l_device = get_component(ElectricLoad, sys, "load1031")

# Decrease the load by increasing the reference value
load_decrease = LoadChange(1.0, l_device, :P_ref, (1.8 * 1.15))

# Define Simulation
sim_dyn_load_decrease = Simulation(
    ResidualModel, #Type of model used
    sys, #system
    pwd(), #folder to output results
    (0.0, 30.0), #time span
    load_decrease, #Type of perturbation
)

#Run the simulation
execute!(
    sim_dyn_load_decrease, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax=0.02, #Maximum step size
)

load_decrease_results = read_results(sim_dyn_load_decrease)

# Plot results
plot_stuff(load_decrease_results, "load_decrease")