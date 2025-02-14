# Set up the package environment
cd(@__DIR__)      # go to the directory of this script
using Pkg         # use the package manager
Pkg.activate(".") # activate the environment defined by the toml files in the parent directory
Pkg.instantiate() # install missing dependencies, make sure environment is ready to use

# Import packages
using PowerSystems

network = System("data/ThreeBusNetwork.raw")

### Buses
# Bus attributes are defined here:
# https://nrel-sienna.github.io/PowerSystems.jl/stable/model_library/generated_ACBus/#ACBus

# Access a single bus by passing the network and the bus number
b1 = get_bus(network, 101)

# Access every bus
for b in get_components(ACBus, network)
    # Set attributes:
    set_voltage_limits!(b, (0.9, 1.1))

    # Query attributes:
    mag = get_magnitude(b)
    print(mag)
    ang = get_angle(b)
    print(ang)
end

### Lines
# Line attributes are defined here:
# https://nrel-sienna.github.io/PowerSystems.jl/stable/model_library/generated_Line/#Line

# Access every line
for l in get_components(Line, network)
    # Set attributes
    set_r!(l, 0)
    set_x!(l, 0.1)

    # Query attributes:
    name = get_name(l)
    println(name)
    # Branch resistance
    r = get_r(l)
    println(r)
    # Branch reactance
    x = get_x(l)
    println(x)
    # Shunt admittance
    b = get_b(l)
    println(b)
    # Angle Limits
    ang_limits = get_angle_limits(l)
    println(ang_limits)
end
