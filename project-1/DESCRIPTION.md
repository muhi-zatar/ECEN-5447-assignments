# Synchronous Machine Model Equations

This document explains the mathematical equations used in the synchronous machine model and their implementation in Julia code. The model is based on the Sauer-Pai approach for power system dynamics and includes a 6th order machine model, IEEE Type-1 AVR, steam turbine with reheat, and governor.

## Table of Contents
- [Reference Frame](#reference-frame)
- [Mechanical Dynamics (Swing Equation)](#mechanical-dynamics-swing-equation)
- [Electrical Dynamics](#electrical-dynamics)
- [Stator Algebraic Equations](#stator-algebraic-equations)
- [IEEE Type-1 AVR Equations](#ieee-type-1-avr-equations)
- [Turbine and Governor Equations](#turbine-and-governor-equations)
- [Network Interface](#network-interface)
- [Initialization Equations](#initialization-equations)

## Reference Frame

The model uses the synchronously rotating d-q reference frame, where:
- The d-axis (direct axis) is aligned with the rotor field
- The q-axis (quadrature axis) is 90 electrical degrees ahead of the d-axis
- The rotor angle δ is the angle between the q-axis and the infinite bus voltage reference

## Mechanical Dynamics (Swing Equation)

### Mathematical Equations

The swing equation describes the rotor dynamics based on the balance of mechanical and electrical torques:

$$\frac{d\delta}{dt} = \omega - \omega_s$$

$$\frac{d\omega}{dt} = \frac{1}{2H}(P_m - P_e - D(\omega - \omega_s))$$

Where:
- δ is the rotor angle in radians
- ω is the rotor speed in per unit
- ω_s is the synchronous speed (1.0 p.u.)
- H is the inertia constant in seconds
- P_m is the mechanical power input
- P_e is the electrical power output
- D is the damping coefficient

### Julia Implementation

```julia
# In the dynamics! function:
du[1] = ω - 1.0  # dδ/dt
du[2] = (Pm - Pe - machine.D*(ω - 1.0))/(2.0*machine.H)  # dω/dt
```

The mechanical power `Pm` is the sum of high-pressure and low-pressure turbine outputs:
```julia
Pm = Pm_hp + Pm_lp
```

The electrical power `Pe` is calculated from terminal voltages and currents:
```julia
Pe = Vd*Id + Vq*Iq + machine.Ra*(Id^2 + Iq^2)
```

## Electrical Dynamics

### Mathematical Equations

The 6th order model includes transient and subtransient dynamics in both d and q axes:

**Transient Dynamics:**

$$\frac{dE'_q}{dt} = \frac{1}{T'_{d0}}(-E'_q - (X_d - X'_d)I_d + V_f)$$

$$\frac{dE'_d}{dt} = \frac{1}{T'_{q0}}(-E'_d + (X_q - X'_q)I_q)$$

**Subtransient Dynamics:**

$$\frac{dE''_q}{dt} = \frac{1}{T''_{d0}}(-E''_q + E'_q + (X'_d - X''_d)I_d)$$

$$\frac{dE''_d}{dt} = \frac{1}{T''_{q0}}(-E''_d + E'_d - (X'_q - X''_q)I_q)$$

Where:
- E'_q, E'_d are the q-axis and d-axis transient voltages
- E''_q, E''_d are the q-axis and d-axis subtransient voltages
- T'_d0, T'_q0 are the open-circuit transient time constants
- T''_d0, T''_q0 are the open-circuit subtransient time constants
- X_d, X_q are the synchronous reactances
- X'_d, X'_q are the transient reactances
- X''_d, X''_q are the subtransient reactances
- I_d, I_q are the stator currents
- V_f is the field voltage

### Julia Implementation

```julia
# Transient dynamics
du[3] = (-Eq_p - (machine.Xd - machine.Xdp) * Id + Vf) / machine.Td0p    # dE'q/dt
du[4] = (-Ed_p + (machine.Xq - machine.Xqp) * Iq) / machine.Tq0p         # dE'd/dt

# Subtransient dynamics
du[5] = (-Eq_pp + Eq_p + (machine.Xdp - machine.Xdpp) * Id) / machine.Td0pp  # dE''q/dt
du[6] = (-Ed_pp + Ed_p - (machine.Xqp - machine.Xqpp) * Iq) / machine.Tq0pp  # dE''d/dt
```

## Stator Algebraic Equations

### Mathematical Equations

The stator voltage equations in the d-q reference frame are:

$$V_d = -R_a I_d + X''_q I_q - E''_d$$

$$V_q = -R_a I_q - X''_d I_d + E''_q$$

These can be rearranged to solve for currents:

$$R_a I_d - X''_q I_q = E''_q - V_d$$

$$X''_d I_d + R_a I_q = E''_q - V_q$$

Which forms a 2×2 system of equations.

The terminal voltage magnitude is:

$$V_t = \sqrt{V_d^2 + V_q^2}$$

### Julia Implementation

```julia
# Matrix solution for currents
A = [machine.Ra -machine.Xqpp; machine.Xdpp machine.Ra]
b = [Eq_pp - Vd; Ed_pp + Vq]
currents = A \ b
Id = currents[1]
Iq = currents[2]

# Terminal voltage magnitude
Vt = sqrt(Vd^2 + Vq^2)
```

## IEEE Type-1 AVR Equations

### Mathematical Equations

The IEEE Type-1 AVR model consists of the following components:

**Voltage Regulator:**

$$\frac{dV_r}{dt} = \frac{1}{T_a}(-V_r + K_a(V_{ref} - V_t - R_f))$$

With limits: $V_{r,min} \leq V_r \leq V_{r,max}$

**Exciter:**

$$\frac{dV_f}{dt} = \frac{1}{T_e}(-K_e V_f + V_r)$$

**Stabilizing Feedback:**

$$\frac{dR_f}{dt} = \frac{1}{T_f}(-R_f + \frac{K_f}{T_f} \frac{dV_f}{dt})$$

Where:
- V_r is the regulator output voltage
- V_f is the field voltage
- R_f is the rate feedback signal
- V_ref is the reference voltage
- V_t is the terminal voltage
- K_a is the regulator gain
- T_a is the regulator time constant
- K_e is the exciter constant
- T_e is the exciter time constant
- K_f is the stabilizer gain
- T_f is the stabilizer time constant

### Julia Implementation

```julia
# Voltage regulator
du[7] = (-Vr + avr.Ka * (Vref - Vt - Rf)) / avr.Ta  # dVr/dt

# Apply limits
if Vr > avr.Vr_max && du[7] > 0
    du[7] = 0.0
elseif Vr < avr.Vr_min && du[7] < 0
    du[7] = 0.0
end

# Exciter
du[8] = (-Vf * avr.Ke + Vr) / avr.Te  # dVf/dt

# Rate feedback
du[9] = (-Rf + avr.Kf * du[8]) / avr.Tf  # dRf/dt
```

## Turbine and Governor Equations

### Mathematical Equations

**Governor Control:**

$$\frac{dP_g}{dt} = \frac{1}{T_g}(-P_g + \frac{\omega_{ref} - \omega}{R})$$

With limits: $P_{min} \leq P_g \leq P_{max}$

**High-Pressure Turbine:**

$$\frac{dP_{m,hp}}{dt} = \frac{1}{T_{ch}}(-P_{m,hp} + F_{hp} P_g)$$

**Low-Pressure Turbine with Reheat:**

$$\frac{dP_{m,lp}}{dt} = \frac{1}{T_{ch}}(-P_{m,lp} + F_{lp} P_g + \frac{P_{m,hp} - F_{hp} P_g}{T_{rh}})$$

Where:
- P_g is the governor output power
- P_m,hp is the high-pressure turbine mechanical power
- P_m,lp is the low-pressure turbine mechanical power
- ω_ref is the speed reference (1.0 p.u. normally)
- ω is the actual rotor speed
- R is the droop coefficient
- T_g is the governor time constant
- T_ch is the steam chest time constant
- T_rh is the reheat time constant
- F_hp, F_lp are the high and low pressure power fractions

### Julia Implementation

```julia
# Governor
du[10] = (-Pg + (ω_ref - ω) / turbine_gov.R) / turbine_gov.Tg  # dPg/dt

# Apply limits
if Pg > turbine_gov.P_max && du[10] > 0
    du[10] = 0.0
elseif Pg < turbine_gov.P_min && du[10] < 0
    du[10] = 0.0
end

# Turbine with reheat
du[11] = (-Pm_hp + turbine_gov.F_hp * Pg) / turbine_gov.T_ch  # dPm_hp/dt
du[12] = (-Pm_lp + turbine_gov.F_lp * Pg + 
          (Pm_hp - turbine_gov.F_hp * Pg) / turbine_gov.T_rh) / turbine_gov.T_ch  # dPm_lp/dt
```

## Network Interface

### Mathematical Equations

The machine connects to the infinite bus through an external impedance. The terminal voltages in the machine reference frame are:

$$V_d = -V_{\infty} \sin(\delta)$$

$$V_q = V_{\infty} \cos(\delta)$$

Where:
- V_∞ is the infinite bus voltage magnitude
- δ is the rotor angle

### Julia Implementation

```julia
# Terminal voltages
Vd = -network.V_∞ * sin(δ)
Vq = network.V_∞ * cos(δ)
```

## Initialization Equations

### Mathematical Equations

For proper initialization, we need to find values for all state variables that satisfy the steady-state condition (all derivatives equal to zero).

**Power Flow Initialization:**

Starting with desired operating point (P_0, Q_0, V_t0):

$$\delta_0 = \sin^{-1}\left(\frac{P_0 X_e}{V_{t0} V_{\infty}}\right)$$

**Initial Currents:**

For a complex power S_0 = P_0 + jQ_0:

$$I_0 = \left(\frac{S_0}{V_{t0}}\right)^*$$

$$I_{d0} = -|I_0| \sin(\theta_0 + \delta_0)$$

$$I_{q0} = |I_0| \cos(\theta_0 + \delta_0)$$

Where θ_0 is the power factor angle.

**Initial Internal Voltages:**

$$E''_{q0} = V_{d0} + R_a I_{d0} + X''_q I_{q0}$$

$$E''_{d0} = V_{q0} - R_a I_{q0} - X''_d I_{d0}$$

$$E'_{q0} = E''_{q0} + (X'_d - X''_d) I_{d0}$$

$$E'_{d0} = E''_{d0} - (X'_q - X''_q) I_{q0}$$

**Field Voltage:**

$$V_{f0} = E'_{q0} + (X_d - X'_d) I_{d0}$$

**AVR Variables:**

$$V_{r0} = V_{f0} K_e$$

$$R_{f0} = \frac{K_f V_{f0}}{T_e}$$

**Governor Variables:**

$$P_{g0} = P_0$$

$$P_{m,hp0} = F_{hp} P_0$$

$$P_{m,lp0} = F_{lp} P_0$$

### Julia Implementation

```julia
# In the initialize! function:

# Calculate initial angle
δ0 = asin(P0 * network.params.X_e / (Vt0 * network.params.V_∞))

# Terminal voltages
Vd0 = -network.params.V_∞ * sin(δ0)
Vq0 = network.params.V_∞ * cos(δ0)

# Calculate currents from desired power
S0 = P0 + im*Q0
I0 = conj(S0/Vt0)
θ0 = angle(I0)

Id0 = -abs(I0) * sin(θ0 + δ0)
Iq0 = abs(I0) * cos(θ0 + δ0)

# Calculate subtransient EMFs
Eq_pp0 = Vd0 + machine.params.Ra * Id0 + machine.params.Xqpp * Iq0
Ed_pp0 = Vq0 - machine.params.Ra * Iq0 - machine.params.Xdpp * Id0

# Calculate transient EMFs
Eq_p0 = Eq_pp0 + (machine.params.Xdp - machine.params.Xdpp) * Id0
Ed_p0 = Ed_pp0 - (machine.params.Xqp - machine.params.Xqpp) * Iq0

# Calculate required field voltage
Vf0 = Eq_p0 + (machine.params.Xd - machine.params.Xdp) * Id0

# AVR equilibrium
Vr0 = Vf0 * avr.params.Ke
Rf0 = avr.params.Kf * Vf0 / avr.params.Te

# Governor equilibrium
Pg0 = P0
Pm_hp0 = turbine_gov.params.F_hp * P0
Pm_lp0 = turbine_gov.params.F_lp * P0
```