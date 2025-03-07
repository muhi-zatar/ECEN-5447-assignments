using DifferentialEquations
using Plots

# Define parameters
g = 9.8                # Gravitational acceleration (m/s^2)
r = 1.0                # Length or scaling factor (m)
ω = 20.0                # Angular velocity (rad/s)
b = 5.0                # Damping coefficient
m = 1.0                # Mass (kg)
ε = (m^2 * g * r) / b^2  # Scaling parameter
T =b/(m*g)
@show r*ω^2/g
# Define the system of ODEs
function system!(du, u, p, t)
    φ, Ω = u
    du[1] = Ω
    du[2] = (1 / ε) * (-sin(φ) + (r * ω^2 / g) * sin(φ) * cos(φ) - Ω)
end

# Initial conditions and time span
u0 = [-π/2, 0.0]        # Initial state [φ, Ω]
tspan = (0.0, 10.0)    # Time span for the simulation

# Solve the ODE
prob = ODEProblem(system!, u0, tspan)
sol = solve(prob, Tsit5())

# Plot the solution
plot(sol, vars=(1), xlabel="Time (s)", ylabel="Variables", label=["φ (Angle)"])
plot!(sol, vars=(2), xlabel="Time (s)", ylabel="Variables", label=["Ω (Angular Velocity)"])
