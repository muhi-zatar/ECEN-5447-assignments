###############################################################################
# three_bus_dynamic.jl
#
# Example code for a dynamic (time-domain) 3-bus, 3-line power system in Julia.
# Each line is modeled as an R-L series with half shunt at each bus. 
# The bus states are the voltages (v1, v2, v3), and the line states 
# are the inductor currents (i12, i13, i23).
#
# This is a single-phase (or "per-phase") model for demonstration.
# You can adapt it to your system's actual data, or to a multi-phase model.
###############################################################################

using DifferentialEquations
using Plots

################################################################################
# 1) Define the ODE function
################################################################################

"""
    three_bus_line_dynamics(du, u, p, t)

Minimal single-phase dynamic model of a 3-bus system with 3 R-L lines and bus+line shunts.

States (u):
  u[1] = v1(t)   [Bus 1 voltage]
  u[2] = v2(t)   [Bus 2 voltage]
  u[3] = v3(t)   [Bus 3 voltage]
  u[4] = i12(t)  [Current on line 1->2]
  u[5] = i13(t)  [Current on line 1->3]
  u[6] = i23(t)  [Current on line 2->3]

Parameters (p) in example:
  p[1] = r12, p[2] = l12, p[3] = b12   (#1 line data)
  p[4] = r13, p[5] = l13, p[6] = b13   (#2 line data)
  p[7] = r23, p[8] = l23, p[9] = b23   (#3 line data)
  p[10] = c1, p[11] = c2, p[12] = c3   (bus self-capacitances, if any)
  p[13] = ig1, p[14] = il1            (net injection/loads at bus1)
  p[15] = ig2, p[16] = il2            (net injection/loads at bus2)
  p[17] = ig3, p[18] = il3            (net injection/loads at bus3)

Units:
  - You can treat r, l, c, b, i, v in consistent SI units (Ohms, Henries, Farads) 
    and base frequency w = 2π·(50 or 60).
  - Or in per-unit. Just be consistent.
"""
dir_raw = "./ThreeBusMultiLoad.raw"
sys = PSY.System(dir_raw)

# Run power flow
pf_result = PF.solve_powerflow(ACPowerFlow{}(), sys)
PSY.set_units_base_system!(sys, "SYSTEM_BASE")

# Get power flow output: complex voltage of each bus terminal 
# (Sienna exports these in pu and rad)
v = pf_result["bus_results"].Vm # [pu-V]
θ = pf_result["bus_results"].θ  # [rad]

# Get power flow output: real (P) and reactive (Q) power injected at each bus
# (Sienna exports these in MW/MVA, so adjust by system base)
P = pf_result["bus_results"].P_net / PSY.get_base_power(sys) # [pu(MW)]
Q = pf_result["bus_results"].Q_net / PSY.get_base_power(sys) # [pu(MVar)]


# -----------------------------------------------------------------------------------------
# Section 2: Calculate relevant quantities needed for running a time domain simulation
# -----------------------------------------------------------------------------------------

##### TODO: Using v,θ,P,Q to find the parameters and initial conditions

# To get you started, here are the branch parameters
PSY.show_components(Line, sys, [:r,:x,:b,:arc]) 

function three_bus_line_dynamics(du, u, p, t)
    # Unpack states
    v1, v2, v3, i12, i13, i23 = u

    # Unpack parameters
    r12, l12, b12,
    r13, l13, b13,
    r23, l23, b23,
    c1, c2, c3,
    ig1, il1,
    ig2, il2,
    ig3, il3 = p

    #---------------------------------------------------------------------------
    # 1) Line currents: 
    #    l12 di12/dt = v1 - v2 - r12*i12
    #    => di12/dt = (v1 - v2 - r12*i12)/l12
    #---------------------------------------------------------------------------
    di12_dt = (v1 - v2 - r12*i12)/l12
    di13_dt = (v1 - v3 - r13*i13)/l13
    di23_dt = (v2 - v3 - r23*i23)/l23

    #---------------------------------------------------------------------------
    # 2) Bus voltages (KCL):
    #
    #    Let "Cbus1 = c1 + b12/2 + b13/2" 
    #       => Cbus1 dv1/dt = ig1 - il1 - i12 - i13
    #    Similarly for bus2 and bus3. 
    #
    #    NOTE: 
    #      i12 flows from bus1->bus2 => sign for bus1 is -i12, 
    #      but sign for bus2 is +i12, etc.
    #---------------------------------------------------------------------------
    Cbus1 = c1 + b12/2 + b13/2
    Cbus2 = c2 + b12/2 + b23/2
    Cbus3 = c3 + b13/2 + b23/2

    dv1_dt = (1/Cbus1)*(ig1 - il1 - i12 - i13)
    dv2_dt = (1/Cbus2)*(ig2 - il2 + i12 - i23)
    dv3_dt = (1/Cbus3)*(ig3 - il3 + i13 + i23)

    # Write results back into du
    du[1] = dv1_dt
    du[2] = dv2_dt
    du[3] = dv3_dt
    du[4] = di12_dt
    du[5] = di13_dt
    du[6] = di23_dt
end

################################################################################
# 2) Choose parameter values for a demonstration
################################################################################

# Example numeric parameters (ARBITRARY for demonstration).
# You should REPLACE these with the actual r,x,b from your .pdf diagram, 
# converting x -> l = x/ω and b -> c = b/(2ω) if you're doing seconds-based simulation.
# Or keep them in per-unit with the assumption that l12 = x12, etc.
r12  = 0.01
l12  = 0.1
b12  = 0.0002

r13  = 0.02
l13  = 0.1
b13  = 0.0001

r23  = 0.015
l23  = 0.08
b23  = 0.00015

# Bus self capacitances (or any net bus-level capacitor). 
# Often, we set c1=c2=c3=0 if we only rely on line-charging. 
# We'll give them small nonzero values for demonstration:
c1   = 0.0005
c2   = 0.0005
c3   = 0.0005

# Net generator currents (ig) and load currents (il) at each bus.
# Could be 0 for no generation or no load. 
ig1, il1 = 1.0, 0.3
ig2, il2 = 1.0, 0.2
ig3, il3 = 0.8, 0.5

# Put them all in a single vector p
p = [
    r12, l12, b12,
    r13, l13, b13,
    r23, l23, b23,
    c1, c2, c3,
    ig1, il1,
    ig2, il2,
    ig3, il3
]

################################################################################
# 3) Initial conditions
################################################################################

# We'll start with all bus voltages at 1.0 p.u. 
# and zero line currents:
u0 = [1.0,  # v1
      1.0,  # v2
      1.0,  # v3
      0.0,  # i12
      0.0,  # i13
      0.0]  # i23

################################################################################
# 4) Build and solve the ODEProblem
################################################################################

tspan = (0.0, 0.2)

# Create the ODEFunction
f = ODEFunction(three_bus_line_dynamics)

# Create the ODEProblem
prob = ODEProblem(f, u0, tspan, p)

# Solve using a solver of your choice, e.g. Tsit5() or Rodas5()
sol = solve(prob, Tsit5())

################################################################################
# 5) Plot results
################################################################################

# Plot all states in one figure
plot(sol, title="3-Bus Dynamic Simulation: All States",
     xlabel="Time [s]")

# Alternatively, separate plots for bus voltages vs line currents
plt1 = plot(sol, idxs=(0, [1,2,3]), 
            label=["v1" "v2" "v3"],
            title="Bus Voltages",
            xlabel="Time [s]")
plt2 = plot(sol, idxs=(0, [4,5,6]), 
            label=["i12" "i13" "i23"],
            title="Line Currents",
            xlabel="Time [s]")
plot(plt1, plt2, layout=(1,2), size=(900,400))
