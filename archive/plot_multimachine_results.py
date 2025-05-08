# Importing libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

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

# Define paths to results
base_path = "../results/Multimachine_Results"

perturbations_to_include = ["line_trip", "load_decrease", "load_increase"]

variations_to_include = ["algebraic_stator_fluxes", "kpq_1.5_rp_0.03"]

# Build dictionary of results
results = {perturbation: {} for perturbation in perturbations_to_include}

for perturbation in perturbations_to_include:
    # Get the baseline results
    base_results_path = f"{base_path}/{perturbation}"
    if os.path.isdir(base_results_path):
        results[perturbation]["baseline"] = pd.read_csv(
            f"{base_results_path}/results.csv"
        )

    # Get any variations
    for variation in variations_to_include:
        variation_path = f"{base_results_path}/{variation}"
        if os.path.isdir(variation_path):
            results[perturbation][variation] = pd.read_csv(
                f"{variation_path}/results.csv"
            )

# --- Angle Plots ---
fig = plt.figure(figsize=(10, 15))

# --- Subplot 1: Load Decrease ---
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(
    results["load_decrease"]["baseline"]["Time"],
    results["load_decrease"]["baseline"]["theta_olc"],
    label="OLC",
)
for key, value in results["load_decrease"].items():
    if key != "baseline":
        ax1.plot(value["Time"], value["theta_olc"], label=key)

ax1.set_title("Load Decrease")
ax1.set_ylabel("Theta (degrees)")
ax1.grid()

# --- Subplot 2 ---
ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(
    results["load_increase"]["baseline"]["Time"],
    results["load_increase"]["baseline"]["theta_olc"],
    label="OLC",
)
for key, value in results["load_increase"].items():
    if key != "baseline":
        ax2.plot(value["Time"], value["theta_olc"], label=key)

ax2.set_title("Load Increase")
ax2.set_ylabel("Theta (degrees)")
ax2.grid()

# --- Subplot 3 ---
ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(
    results["line_trip"]["baseline"]["Time"],
    results["line_trip"]["baseline"]["theta_olc"],
    label="OLC",
)
for key, value in results["line_trip"].items():
    if key != "baseline":
        ax3.plot(value["Time"], value["theta_olc"], label=key)

ax3.set_title("Line Trip")
ax3.set_ylabel("Theta (degrees)")
ax3.grid()

ax3.set_xlabel("Time (s)")

# --- Legend ---
lines, labels = [], []
for ax in [ax1]:
    line, label = ax.get_legend_handles_labels()
    lines += line
    labels += label

fig.legend(lines, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3)
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig("angles.png", bbox_inches="tight")
