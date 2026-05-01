"""Generates a density plot for side-chain spatial angles in XRX helical motifs,
categorized by the steric bulk of the preceding amino acid.

Replaces R-based plotting with matplotlib for seamless pipeline integration.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib

# Force matplotlib to not use any Xwindows backend
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Define categories and their corresponding colors
VOLUME_TIERS = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]
PALETTE = {
    "Tiny": "#e6e2d3",
    "Small": "#f2c879",
    "Intermediate": "#f28e5c",
    "Large": "#e4572e",
    "Bulky": "#b11226",
}


def generate_angle_plot():
    # Fetch file paths from arguments or use defaults
    input_tsv = sys.argv[1] if len(sys.argv) > 1 else "final/angles.tsv"
    output_image = sys.argv[2] if len(sys.argv) > 2 else "final/ss_profile_HHH_for_arg_with_valid_runs.png"

    # Load data
    angle_data = pd.read_csv(input_tsv, sep="\t")

    # Normalize angles to strictly fall within the [-180, 180] degree range
    angle_data["angle"] = ((angle_data["angle"] + 180) % 360) - 180
    
    # Filter to keep only the recognized size classifications
    angle_data = angle_data[angle_data["size_class"].isin(VOLUME_TIERS)]

    total_samples = len(angle_data)

    # Initialize the plot canvas with the designated gray background
    figure, axes = plt.subplots(figsize=(10, 6))
    figure.patch.set_facecolor("#bdbdbd")
    axes.set_facecolor("#bdbdbd")

    # Create 1000 evenly spaced points between -180 and 180 for smooth KDE curves
    x_axis_range = np.linspace(-180, 180, 1000)

    # Plot the Kernel Density Estimate (KDE) for each volume tier
    for tier in VOLUME_TIERS:
        tier_angles = angle_data[angle_data["size_class"] == tier]["angle"].values
        
        # Need at least 2 data points to compute a density curve
        if len(tier_angles) < 2:
            continue
            
        density_estimator = gaussian_kde(tier_angles, bw_method="silverman")
        y_density = density_estimator(x_axis_range)
        axes.plot(x_axis_range, y_density, color=PALETTE[tier], linewidth=1.6, label=tier)

        # Configure axes ranges and tick marks
    axes.set_xlim(-180, 180)
    axes.set_xticks(np.arange(-180, 181, 50))

    # Overlay vertical dashed blue grid lines
    for tick_val in np.arange(-180, 181, 50):
        axes.axvline(tick_val, color="blue", linestyle=":", linewidth=0.5, alpha=0.7)

    # Apply axis labels and plot title
    axes.set_xlabel(r"Angle between adjacent C-$\alpha$ $\to$ Centroid vectors [°]")
    axes.set_ylabel("Norm. Freq. [A.U.]")
    axes.set_title(f"Tripeptide (XRX) in Helix (n = {total_samples})")

    # Format the legend box
    legend_box = axes.legend(loc="upper left", facecolor="white", edgecolor="black")
    legend_box.get_frame().set_linewidth(0.8)

    # Clean up the visual borders by hiding the top and right spines
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)

    # Ensure the target directory exists and save the figure
    target_dir = os.path.dirname(output_image) or "."
    os.makedirs(target_dir, exist_ok=True)
    
    plt.tight_layout()
    plt.savefig(output_image, dpi=300, facecolor=figure.get_facecolor())
    
    print(f"Plot successfully saved to: {output_image}")
    print(f"Total data points plotted (n) = {total_samples}")


if __name__ == "__main__":
    generate_angle_plot()