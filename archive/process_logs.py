# Importing libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set font size
plt.rcParams.update(
    {
        "axes.titlesize": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 12.5,
        "ytick.labelsize": 12.5,
        "legend.fontsize": 16,
    }
)

# Read the logs
logs_line = pd.read_csv("inverter_simulation_results_line_trip.csv")
logs_115 = pd.read_csv("inverter_simulation_results_load_115.csv")
logs_85 = pd.read_csv("inverter_simulation_results_load_85.csv")

# Read the time
time_line = logs_line["Time"]
time_115 = logs_115["Time"]
time_85 = logs_85["Time"]

# --- Theta PLL and OLC ---
theta_pll_line = logs_line["theta_pll"]
theta_pll_115 = logs_115["theta_pll"]
theta_pll_85 = logs_85["theta_pll"]
theta_olc_line = logs_line["theta_olc"]
theta_olc_115 = logs_115["theta_olc"]
theta_olc_85 = logs_85["theta_olc"]

fig = plt.figure(figsize=(10, 15))
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(time_115, theta_pll_115, label="PLL", color="tab:brown")
ax1.plot(time_115, theta_olc_115, label="OLC", color="tab:cyan")
ax1.set_title("Load Decrease")
ax1.set_ylabel("Theta (degrees)")
ax1.grid()

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(time_85, theta_pll_85, label="PLL", color="tab:brown")
ax2.plot(time_85, theta_olc_85, label="OLC", color="tab:cyan")
ax2.set_title("Load Increase")
ax2.set_ylabel("Theta (degrees)")
ax2.grid()

ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(time_line, theta_pll_line, label="PLL", color="tab:brown")
ax3.plot(time_line, theta_olc_line, label="OLC", color="tab:cyan")
ax3.set_title("Line Trip")
ax3.set_ylabel("Theta (degrees)")
ax3.set_xlabel("Time (s)")
ax3.grid()

handles, labels = ax3.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2)
fig.tight_layout(rect=[0, 0.025, 1, 1])
fig.savefig("theta_pll_olc.png", bbox_inches="tight")

# --- Voltage Magnitudes ---
v_mag_line_1 = logs_line["voltage_magnitude_1"]
v_mag_line_2 = logs_line["voltage_magnitude_2"]
v_mag_line_3 = logs_line["voltage_magnitude_3"]
v_mag_115_1 = logs_115["voltage_magnitude_1"]
v_mag_115_2 = logs_115["voltage_magnitude_2"]
v_mag_115_3 = logs_115["voltage_magnitude_3"]
v_mag_85_1 = logs_85["voltage_magnitude_1"]
v_mag_85_2 = logs_85["voltage_magnitude_2"]
v_mag_85_3 = logs_85["voltage_magnitude_3"]

fig = plt.figure(figsize=(10, 15))
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(time_115, v_mag_115_1, label="Infinite Bus", color="tab:blue")
ax1.plot(time_115, v_mag_115_2, label="Converter", color="tab:orange")
ax1.plot(time_115, v_mag_115_3, label="Load", color="magenta")
ax1.set_title("Load Decrease")
ax1.set_ylabel("Voltage Magnitude (p.u.)")
ax1.grid()

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(time_85, v_mag_85_1, label="Infinite Bus", color="tab:blue")
ax2.plot(time_85, v_mag_85_2, label="Converter", color="tab:orange")
ax2.plot(time_85, v_mag_85_3, label="Load", color="magenta")
ax2.set_title("Load Increase")
ax2.set_ylabel("Voltage Magnitude (p.u.)")
ax2.grid()

ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(time_line, v_mag_line_1, label="Infinite Bus", color="tab:blue")
ax3.plot(time_line, v_mag_line_2, label="Converter", color="tab:orange")
ax3.plot(time_line, v_mag_line_3, label="Load", color="magenta")
ax3.set_title("Line Trip")
ax3.set_ylabel("Voltage Magnitude (p.u.)")
ax3.set_xlabel("Time (s)")
ax3.grid()

handles, labels = ax3.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3)
fig.tight_layout(rect=[0, 0.025, 1, 1])
fig.savefig("voltage_magnitude.png", bbox_inches="tight")


# --- Line Currents ---
i_line_12 = np.sqrt(logs_line["I_12_r"] ** 2 + logs_line["I_12_i"] ** 2)
i_line_13 = np.sqrt(logs_line["I_13_r"] ** 2 + logs_line["I_13_i"] ** 2)
i_line_23 = np.sqrt(logs_line["I_23_r"] ** 2 + logs_line["I_23_i"] ** 2)
i_115_12 = np.sqrt(logs_115["I_12_r"] ** 2 + logs_115["I_12_i"] ** 2)
i_115_13 = np.sqrt(logs_115["I_13_r"] ** 2 + logs_115["I_13_i"] ** 2)
i_115_23 = np.sqrt(logs_115["I_23_r"] ** 2 + logs_115["I_23_i"] ** 2)
i_85_12 = np.sqrt(logs_85["I_12_r"] ** 2 + logs_85["I_12_i"] ** 2)
i_85_13 = np.sqrt(logs_85["I_13_r"] ** 2 + logs_85["I_13_i"] ** 2)
i_85_23 = np.sqrt(logs_85["I_23_r"] ** 2 + logs_85["I_23_i"] ** 2)

fig = plt.figure(figsize=(10, 15))
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(time_115, i_115_12, label="I_12", color="blue")
ax1.plot(time_115, i_115_13, label="I_13", color="tab:green")
ax1.plot(time_115, i_115_23, label="I_23", color="tab:purple")
ax1.set_title("Load Decrease")
ax1.set_ylabel("Line Current (p.u.)")
ax1.grid()

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(time_85, i_85_12, label="I_12", color="blue")
ax2.plot(time_85, i_85_13, label="I_13", color="tab:green")
ax2.plot(time_85, i_85_23, label="I_23", color="tab:purple")
ax2.set_title("Load Increase")
ax2.set_ylabel("Line Current (p.u.)")
ax2.grid()

ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(time_line, i_line_12, label="I_12", color="blue")
ax3.plot(time_line, i_line_13, label="I_13", color="tab:green")
ax3.plot(time_line, i_line_23, label="I_23", color="tab:purple")
ax3.set_title("Line Trip")
ax3.set_xlabel("Time (s)")
ax3.set_ylabel("Line Current (p.u.)")
ax3.grid()

handles, labels = ax3.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3)
fig.tight_layout(rect=[0, 0.025, 1, 1])
fig.savefig("line_current.png", bbox_inches="tight")


# --- Active and Reactive Power ---
p_line_1, p_line_2, p_line_3 = logs_line["P_1"], logs_line["P_2"], logs_line["P_3"]
q_line_1, q_line_2, q_line_3 = logs_line["Q_1"], logs_line["Q_2"], logs_line["Q_3"]
p_115_1, p_115_2, p_115_3 = logs_115["P_1"], logs_115["P_2"], logs_115["P_3"]
q_115_1, q_115_2, q_115_3 = logs_115["Q_1"], logs_115["Q_2"], logs_115["Q_3"]
p_85_1, p_85_2, p_85_3 = logs_85["P_1"], logs_85["P_2"], logs_85["P_3"]
q_85_1, q_85_2, q_85_3 = logs_85["Q_1"], logs_85["Q_2"], logs_85["Q_3"]

fig = plt.figure(figsize=(10, 15))
ax1 = fig.add_subplot(3, 1, 1)
(l1,) = ax1.plot(time_115, p_115_1, label="Infinite Bus P", color="tab:blue")
(l2,) = ax1.plot(time_115, p_115_2, label="Converter P", color="tab:orange")
(l3,) = ax1.plot(time_115, p_115_3, label="Load P", color="magenta")
(l4,) = ax1.plot(
    time_115, q_115_1, label="Infinite Bus Q", linestyle="--", color="tab:blue"
)
(l5,) = ax1.plot(
    time_115, q_115_2, label="Converter Q", linestyle="--", color="tab:orange"
)
(l6,) = ax1.plot(time_115, q_115_3, label="Load Q", linestyle="--", color="magenta")
ax1.set_title("Load Decrease")
ax1.set_ylabel("Power (p.u.)")
ax1.grid()

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(time_85, p_85_1, color="tab:blue")
ax2.plot(time_85, p_85_2, color="tab:orange")
ax2.plot(time_85, p_85_3, color="magenta")
ax2.plot(time_85, q_85_1, linestyle="--", color="tab:blue")
ax2.plot(time_85, q_85_2, linestyle="--", color="tab:orange")
ax2.plot(time_85, q_85_3, linestyle="--", color="magenta")
ax2.set_title("Load Increase")
ax2.set_ylabel("Power (p.u.)")
ax2.grid()

ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(time_line, p_line_1, color="tab:blue")
ax3.plot(time_line, p_line_2, color="tab:orange")
ax3.plot(time_line, p_line_3, color="magenta")
ax3.plot(time_line, q_line_1, linestyle="--", color="tab:blue")
ax3.plot(time_line, q_line_2, linestyle="--", color="tab:orange")
ax3.plot(time_line, q_line_3, linestyle="--", color="magenta")
ax3.set_title("Line Trip")
ax3.set_xlabel("Time (s)")
ax3.set_ylabel("Power (p.u.)")
ax3.grid()

# Grouped legend handles and labels
handles = [l1, l4, l2, l5, l3, l6]
labels = [
    "Infinite Bus P",
    "Infinite Bus Q",
    "Converter P",
    "Converter Q",
    "Load P",
    "Load Q",
]

fig.legend(
    handles,
    labels,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.0),
    ncol=3,
    columnspacing=1.2,
)
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig("power.png", bbox_inches="tight")
