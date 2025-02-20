cd(@__DIR__)
using Pkg
Pkg.activate(".")
#Pkg.resolve()
#Pkg.instantiate()

# Import packages
using DifferentialEquations
using Plots

### Based on Example 1.9 from Milano's Eigenvalue Problems in Power Systems

# Define function to be used by integration method  
function omib_classical!(du, u, p, t)
    # Rename states for clarity
    δr = u[1]  # [rad] rotor angular position 
    ωr = u[2]  # [pu(rad/s)] rotor angular speed 

    # Rename parameters for clarity
    M1 = p[1]   # [s*pu(MW)]
    D1 = p[2]   # [pu(MW)]
    erq1 = p[3] # [pu(kV)]
    Pm1 = p[4]  # [pu(MW)]
    ω0 = p[5]   # [rad/s]
    vth = p[6]  # [pu(kV)]
    θth = p[7]  # [rad]
    Xtot = p[8] # [pu(Ω)]

    # Define the differential equations
    #  u[1] = δr, so du[1] = d(δr)/dt
    #  u[2] = ωr, so du[2] = d(ωr)/dt
    du[1] = ω0 * (ωr - 1)
    du[2] = (1/M1) * ( Pm1 - ((erq1*vth/Xtot)*sin(δr - θth)) - (D1*(ωr - 1)) )

end

# Define parameters
M1 = 234    # [s*pu(MW)] mechanical starting time 
D1 = 0      # [pu(MW)] damping coefficient of 
erq1 = 1.05 # [pu(kV)] EMF of generator behind the reactance
Pm1 = 6.67  # [pu(MW)] mechanical power 
ω0 = 377    # [rad/s] reference angular frequency (speed of grid reference frame)
vth = 1.0   # [pu(kV)] voltage of infinite bus
θth = 0     # [rad] phase angle of infinite bus
Xtot = 0.1  # [pu(Ω)] internal X of mach. + line X from mach. to IB + Thev. X of grid
p = [M1, D1, erq1, Pm1, ω0, vth, θth, Xtot]

# Find fixed points (which we can conveniently calculate for this particular system)
δr_fp_stable = asin(Pm1 * Xtot/(erq1*vth)) - θth 
δr_fp_unstable = pi - asin(Pm1 * Xtot/(erq1*vth)) - θth 
ωr_fp = 1.0

# # Check that the fixed points are valid
du0 = Vector{Float64}(undef,2) # return object for derivative calculation
omib_classical!(du0, [δr_fp_stable; ωr_fp], p, 0)
omib_classical!(du0, [δr_fp_unstable; ωr_fp], p, 0)

# Define initial conditions
ϵ = 0.001 # rotor angle perturbation
δr0 = δr_fp_stable + ϵ  # [rad] rotor angular position 
ωr0 = ωr_fp             # [pu(rad/s)] rotor angular speed
u0 = [δr0, ωr0]

# Define length of simulation
tspan = (0.0, 100.0) # [s]

# Build problem
prob = ODEProblem(omib_classical!, u0, tspan, p)

# Run simulation
sol = solve(prob)

# Plot the time-domain trajectory of the state variables on subplots
p = plot(sol, 
    layout = (2,1), 
    label=["δr" "ωr"],
    suptitle="Time Domain: ϵ = $(ϵ)",
    title=["Rotor Angular Position" "Rotor Angular Velocity  "], 
    xlabel= "Time [s]", 
    ylabel = ["Rotor Angular Position [rad]" "Rotor Angular Velocity [pu(rad/s)]"], 
    color=["blue" "red"]
    )
hline!(p[1], [δr_fp_stable], label="Fixed Point", color=:black)
hline!(p[2], [ωr_fp], label="Fixed Point", color=:black)

# Plot the phase portrait of the two state variables
plot(sol, 
    idxs=(1,2),  # 0 = time, 1 = first state var, 2 = second state var, ...
    title="Phase Portrait: Velocity vs Position", 
    xlabel="Rotor Angular Position [rad]", 
    ylabel="Rotor Angular Velocity [pu(rad/s)]"
    )

# Make plot from book to visualize fixed points.
Pe(δr) = (erq1*vth/Xtot)*sin(δr - θth) # external torque [Nm]
δr_range = collect(range(0, pi, length=100))
plot(δr_range, 
    Pe.(δr_range), 
    title="Sychronous Generator P-delta Relationship", 
    xlabel="δr [rad]", 
    ylabel="P [pu(MW)]", 
    label="Pe(δr)"
    )
hline!([Pm1], label="Pm")

# Show example of integration steps vs interpolated values
plot(sol.t, 
    sol[1,:], 
    title="Example of Integration Steps vs. Interpolated Values",
    label="Exact Integration Steps", 
    marker=:circle
    ) 
plot!(sol, 
    idxs=(0,1), 
    label="Interpolated Curve",
    )
