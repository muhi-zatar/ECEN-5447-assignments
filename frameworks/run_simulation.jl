cd(@__DIR__)
using Pkg
Pkg.activate("../.")
#Pkg.resolve() # make sure Manifest matches Project
#Pkg.instantiate() # install missing dependencies

# Importing necessary Modules
using PowerSystemCaseBuilder
using PowerSimulationsDynamics
using PowerNetworkMatrices
using Plots
using Sundials
using InfrastructureSystems
using PowerSystems
using Logging
const PSY = PowerSystems

# Build the system frp, PSID systems, this should be similar to ours.
sys = build_system(PSIDTestSystems, "psid_test_threebus_multimachine_dynlines")

# See which buses have generators
gens = collect(get_components(ThermalStandard, sys))
ref_gen = gens[2]
machine_gen = gens[1]

# Remove dynamic injectors (to prepare for replacement by our model)
remove_component!(sys, get_dynamic_injector(ref_gen))
remove_component!(sys, get_dynamic_injector(machine_gen))

# Remove the ThermalStandard generator that's currently at the reference bus
remove_component!(sys, ref_gen)

### Add an infinite voltage source
# Define the slack bus
slack_bus = first(get_components(x -> get_bustype(x) == ACBusTypes.REF, Bus, sys))

# Define an infinite source
inf_source = Source(;
    name="InfBus", #name
    available=true, #availability
    active_power=0.0,
    reactive_power=0.0,
    bus=slack_bus, #bus
    R_th=0.0, #Rth
    X_th=0.1, #Xth
);

# Add the infinite source
add_component!(sys, inf_source)

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
    # γ_d1 = ,
    # γ_q1 = ,
    # γ_d2 = ,
    # γ_q2 = ,
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
exciter = PSY.EXST1(
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
    avr=exciter,
    prime_mover=governor,
    # exciter = exciter,
    pss=PSY.PSSFixed(V_pss=0.0),
    base_power=get_base_power(machine_gen),
    # InfrastructureSystems.InfrastructureSystemsInternal()
)

# Add our generator to the system
add_component!(sys, dg1, machine_gen)


### Run a Dynamic Simulation ###
# Define Simulation
tspan = (0.0, 30.0)
sim = Simulation(
    ResidualModel,
    sys,
    pwd(),
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

# Plot some things
V_mag = get_voltage_magnitude_series(results, 102)
plot(V_mag)

### Line Trip Perturbations

#Make a copy of the original system
new_sys = deepcopy(sys);

#Remove Line "BUS 1-BUS 2-i_1"
remove_component!(DynamicBranch, new_sys, "BUS 2-BUS 3-i_1")

# Weird update of Ybus – we don't know why the tutorial does this
# fault_branches2 = get_components(DynamicBranch, new_sys)

# for br in fault_branches2
#     if get_name(br) == "BUS 1-BUS 3-i_1"
#         br.r = 3 * br.r
#         br.x = 3 * br.x
#         b_new = (from=br.b.from / 3, t=br.b.to / 3)
#         br.b = b_new
#     end
# end

# Obtain the new Ybus
Ybus_fault_dyn = Ybus(new_sys).data

# Hack – set our Ybus to exactly the same values as in Step 3.1 of https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/tutorials/tutorial_dynamic_lines/
# Ybus_fault_dyn[1, 1] = 0.91954 - im * 10.9011
# Ybus_fault_dyn[2, 2] = 0.689655 - im * 8.17586
# Ybus_fault_dyn[1, 3] = -0.229885 + im * 2.75862
# Ybus_fault_dyn[3, 1] = -0.229885 + im * 2.75862
# Ybus_fault_dyn[3, 3] = 0.229885 - im * 2.72529

Ybus_change_dyn = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault_dyn, #New YBus
)

# Define Simulation
sim_dyn = Simulation(
    ResidualModel, #Type of model used
    sys, #system
    pwd(), #folder to output results
    (0.0, 30.0), #time span
    Ybus_change_dyn, #Type of perturbation
)

#Run the simulation
execute!(
    sim_dyn, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax=0.02, #Maximum step size
)

line_trip_results = read_results(sim_dyn)

# Plot some things
V_mag = get_voltage_magnitude_series(line_trip_results, 102)
plot(V_mag)
