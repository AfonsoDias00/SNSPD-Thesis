import numpy as np
import matplotlib.pyplot as plt
import tdgl
from scipy.interpolate import griddata
from pathlib import Path

# Setup paths based on your existing script structure
BASE_DIR = Path(__file__).resolve().parent
H5_PATH = BASE_DIR / "results_four_holes" / "h5" / "four_symmetric_B0mT.h5"
FIGURE_DIR = BASE_DIR / "results_four_holes" / "figures"

def main():
    if not H5_PATH.exists():
        raise FileNotFoundError(f"Cannot find {H5_PATH}. Run the 4-hole simulation first.")

    # 1. Load the existing solution
    print(f"Loading {H5_PATH.name}...")
    solution = tdgl.Solution.from_hdf5(str(H5_PATH))
    device = solution.device

    # 2. Extract mesh coordinates and the final current density snapshot
    # device.points contains the (x, y) coordinates of every mesh site
    x = device.points[:, 0]
    y = device.points[:, 1]
    
    # solution.current_density returns the (J_x, J_y) vectors at each site
    J_vectors = solution.current_density
    J_y = J_vectors[:, 1]  # We only care about the longitudinal current component

    # 3. Create a high-resolution 1D line across the width of the wire at y = 0
    # The wire width is 0.2 um, so from x = -0.1 to 0.1
    line_x = np.linspace(-0.1, 0.1, 1000)
    line_y = np.zeros_like(line_x)  # y = 0 (cutting exactly through the top holes)

    # 4. Interpolate the unstructured mesh data onto our 1D line
    J_y_interp = griddata((x, y), J_y, (line_x, line_y), method="cubic")

    # 5. Plot the 1D current density profile
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Plot the absolute current density
    ax.plot(line_x, np.abs(J_y_interp), color="tab:blue", linewidth=2, label="|J_y| profile")

    # Visually mark the hole locations (hole diameter is 50 nm, centers at x = -0.055 and +0.055)
    # The holes span from center - 0.025 to center + 0.025
    ax.axvspan(-0.080, -0.030, color="gray", alpha=0.3, label="Hole voids")
    ax.axvspan(0.030, 0.080, color="gray", alpha=0.3)

    ax.set_title("1D Current Density Profile across Symmetric 4-Hole Grid (y=0)")
    ax.set_xlabel("Wire Width: x [μm]")
    ax.set_ylabel("Current Density $|J_y|$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper center")
    
    # Save the 1D profile plot
    out_file = FIGURE_DIR / "four_symmetric_B0mT_1D_profile.png"
    fig.savefig(out_file, bbox_inches="tight")
    print(f"Saved 1D profile to {out_file}")
    

if __name__ == "__main__":
    main()