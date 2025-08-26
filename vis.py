"""Functions for visualizing results in terminal and as plots.

"""

# Standard library
import time

# Third party
import numpy as np
import matplotlib.pyplot as plt

# Local
import calc

def print_ten_first_energies(N: int, L: int) -> None:
    """Calculates and prints first 10 energies.
    
    """
    # Calculate and time numerical calculations
    start = time.perf_counter()
    approx_energies = calc.numerical_energies(N, L)
    end = time.perf_counter()

    exec_time = end - start

    # We only care about the 10 first energies
    actual_energies = calc.analytical_energies(10)
    approx_energies = approx_energies[:10]

    # Find error percentages of the approximated energies
    errors = np.abs(approx_energies/actual_energies - 1)*100

    # Print detailed information
    print()
    print(f"         N = {N} | L = {L} | exec. time: {exec_time*1000:.3f} ms")  
    print()
    print("|  State | Actual Energy | Approximate Energy |   Error %  |")
    print("|--------|---------------|--------------------|------------|")
    for i in range(10):
        print(f"| n = {i + 1:>2} | {actual_energies[i]:>13} "
              f"| {approx_energies[i]:>18.8f} | {errors[i]:>10.5f} |")
    print()

def plot_N_vs_L(N_min: int = 12, N_max: int = 20, 
                L_min: int =  8, L_max: int = 20) -> None:
    """Plots N vs L plot where color indicates all first 10 energies
    
    """
    # Set simulation values
    L_values = np.arange(L_min, L_max + 1)
    N_values = np.arange(N_min, N_max + 1)

    # Create value grid
    L_gird, N_grid = np.meshgrid(L_values, N_values)

    # Calculate mean errors
    error_grid = np.zeros((len(N_values), len(L_values)))
    
    actual_energies = calc.analytical_energies(10)
    for i, N in enumerate(N_values):
        for j, L in enumerate(L_values):
            approx_energies = calc.numerical_energies(N, L)[:10]
            energy_error = np.abs(approx_energies/actual_energies - 1)*100
            
            if np.all(energy_error < 1):
                error_grid[i, j] = np.mean(energy_error)
            else:
                error_grid[i, j] = np.nan  # Not interested in values > 1% error

    # Create custom color map where "white" is used for all np.nan values
    cmap = plt.get_cmap("winter")
    cmap.set_bad(color="white")

    # Draw color plot
    color = plt.pcolormesh(L_gird, N_grid, error_grid, cmap=cmap)
    plt.colorbar(color, label="Mean error %")

    # Set ticks and titles
    plt.xticks(L_values)
    plt.yticks(N_values)
    plt.xlabel("L")
    plt.ylabel("N", rotation=0)
    plt.title("N and L combinations that produce <1% error\n"
              "for first 10 energy states")

    # Draw custom grid lines (plt.grid() does not support offsetting from ticks)
    gird_vlines = L_values - 0.5
    grid_hlines = N_values - 0.5
    plt.gca().vlines(x=gird_vlines, 
                     ymin=np.min(grid_hlines), 
                     ymax=np.max(grid_hlines + 1),
                     color="Gray", linewidth=0.5)
    plt.gca().hlines(y=grid_hlines,
                     xmin=np.min(gird_vlines), 
                     xmax=np.max(gird_vlines + 1), 
                     color="Gray", linewidth=0.5)
    plt.show()