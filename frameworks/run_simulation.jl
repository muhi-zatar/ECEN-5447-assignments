# Importing necessary Modules
using PowerSystemCaseBuilder
using PowerSimulationsDynamics
using PowerNetworkMatrices
using Plots
using Sundials
using InfrastructureSystems
using PowerSystems
const PSY = PowerSystems

# Build the system frp, PSID systems, this should be similar to ours.
sys = build_system(PSIDTestSystems, "psid_test_threebus_multimachine_dynlines")

# See which buses have generators
gens = collect(get_components(ThermalStandard, sys))
gen1 = gens[1]
gen2 = gens[2]


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
    Dict("S10" => 0.0, "S12" => 0.0),
    3.148, # H
    2.0, # D
    1.0,    # ω_ref (typically 1.0 pu)
    60.0,   # frequency (Hz)
    [:δ, :ω, :Eqp, :Edp, :psikd, :psikq],  # state names
    6,      # number of states
    InfrastructureSystems.InfrastructureSystemsInternal()
)

# Define Shaft with single mass
shaft_model = SingleMass(    
    H = 3.148,
    D = 2.0,
    )

# Degine GasTG governer
governor = GasTG(
    0.05,    # R
    0.2,     # T1
    0.4,     # T2
    0.04,    # T3
    1.0,     # AT
    2.0,     # Kt
    (0.0, 6.0), # Vlim
    0.0    # D_turb
)

# Define EXST1 AVR
exciter = PSY.EXST1(
    0.001, # Tr
    (0.0, 6.0), # Vi_lim
    400.0, # Tc
    0.02, # Tb
    200, # Ka
    0.8, # Ta
    (-6.0, 6.0), # Vr_limit
    0.0, # Kc
    0.0, # Kf
    0.0, # Tf
    1.0, # V_ref
)

# Define dynamic generator wrapper
dg1 = DynamicGenerator(
    name = "generator-101-1",
    ω_ref = 1.0,
    machine = sp_machine,
    shaft = shaft_model,
    avr = exciter,
    prime_mover = governor,
    # exciter = exciter,
    pss = PSY.PSSFixed(V_pss = 0.0),
    base_power = get_base_power(gen2),    
    # InfrastructureSystems.InfrastructureSystemsInternal()
)

# Define Simulationtspan = (0.0, 30.0)
try
    sim = Simulation(
        ResidualModel,
        sys,
        pwd(),
        tspan
    )
catch e
    println("Simulation error: ", e)
    println(stacktrace(catch_backtrace()))
end

# execute the simulation
execute!(
    sim, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = 0.02, #Maximum step size
    )

results = read_results(sim)