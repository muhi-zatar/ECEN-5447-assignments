# Set up the package environment
cd(@__DIR__)      # go to the directory of this script
using Pkg         # use the package manager
Pkg.activate(".") # activate the environment defined by the toml files in the parent directory
Pkg.instantiate() # install missing dependencies, make sure environment is ready to use

# Import packages
using PowerSystems

network = System("data/ThreeBusMultiLoad.raw")

# Print a description of the network
print(network)

# Show the different static components
show_components(network, ACBus)
show_components(network, Arc)
show_components(network, StandardLoad)
show_components(network, ThermalStandard)

### Buses
# Bus attributes are defined here:
# https://nrel-sienna.github.io/PowerSystems.jl/stable/model_library/generated_ACBus/#ACBus

# Access a single bus
bus_1 = get_component(ACBus, network, "BUS 1")
bus_2 = get_component(ACBus, network, "BUS 2")
bus_3 = get_component(ACBus, network, "BUS 3")

# Access every bus
for b in get_components(ACBus, network)
    # Set attributes:
    #set_voltage_limits!(b, (0.9, 1.1))

    # Query attributes:
    println("Bus: $(get_name(b))\n")
    println("Bus number: $(get_number(b))\n")
    println("Type = $(get_bustype(b))\n")
    println("Base voltage: $(get_base_voltage(b))\n")
    #mag = get_magnitude(b)
    #print(mag)
    #ang = get_angle(b)
    #print(ang)
end

### Lines
# Line attributes are defined here:
# https://nrel-sienna.github.io/PowerSystems.jl/stable/model_library/generated_Line/#Line

# Access a single line
line_12 = get_component(Line, network, "BUS 1-BUS 2-i_1")
line_13 = get_component(Line, network, "BUS 1-BUS 3-i_1")
line_23 = get_component(Line, network, "BUS 2-BUS 3-i_1")

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

### Loads
# Load attributes are defined here:
#https://nrel-sienna.github.io/PowerSystems.jl/stable/model_library/generated_StandardLoad/

# Access a single load
load_1 = get_component(StandardLoad, network, "load1011")

# Access every load
for ld in get_components(StandardLoad, network)
    # Set attributes
    set_available!(ld, false)
    set_max_constant_active_power!(ld, 2.0)
    set_constant_active_power!(ld, 2.0)

    # Query attributes
    name = get_name(ld)
    println(name)
    bus_name = get_name(get_bus(ld))
    println(bus_name)
    base_power = get_base_power(ld)
    println(base_power)
end


### Generators
# Generator attributes depend on the type of generation. In this case, the network only has type 
# "ThermalStandard" generators. Attributes for ThermalStandard generators are defined here:
# https://nrel-sienna.github.io/PowerSystems.jl/stable/model_library/generated_ThermalStandard/

# Access a single generator
gen_1 = get_component(ThermalStandard, network, "generator-102-1")

# Access every generator
for g in get_components(ThermalStandard, network)
    # Set attributes
    set_available!(g, false)
    set_active_power!(g, 1.1)
    set_active_power_limits!(g, (0.5, 1.5))

    # Query attributes
    base_power = get_base_power(g)
    println(base_power)
    rating = get_rating(g)
    println(rating)
    ramp_limits = get_ramp_limits(g)
    println(ramp_limits)

end


### Per-unit
# Per-unit conventions are described here:
# https://nrel-sienna.github.io/PowerSystems.jl/stable/explanation/per_unit/

# Find the units base for the system
units_base = get_units_base(network)
print(units_base)

# Print the base power in the system units base
system_base_power = get_base_power(network)
println(system_base_power)

