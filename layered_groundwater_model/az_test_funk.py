import numpy as np
import matplotlib as mpl; mpl.use("WebAgg")
import matplotlib.pyplot as plt
import flopy
import pandas as pd

def plot_heads_and_concentrations(scenario=2, show_plots=True):
    """
    Load MODFLOW 6 simulation data and plot head and concentration time series
    for wells defined in a CSV file.

    Parameters:
    - scenario (int): Scenario number to load data from.
    - show_plots (bool): Whether to display the plots.
    """
    file_path = f"./outputs/s{scenario}/"

    # Load simulation
    mf = flopy.mf6.MFSimulation.load(sim_name="mfsim", sim_ws=file_path)
    gwf = mf.get_model(f"gwf_s{scenario}")
    gwt = mf.get_model(f"gwt_s{scenario}")

    # Extract head and concentration data
    hds = gwf.output.head().get_alldata()
    ucn = gwt.output.concentration().get_alldata()

    print("Grabbed head and concentration data.")

    # Load well information
    well = pd.read_csv('./inputs/well_info.csv')
    coords = well.set_index('type')
    x_bin = coords.loc['well']['x_coord'] // gwf.modelgrid.delr[0]
    z_bin = coords.loc['well']['z_coord'] // gwf.modelgrid.delz[0, 0, 0]

    # Plot heads
    plt.figure(1)
    for z, x, name in zip(z_bin, x_bin, coords.loc['well']['id']):
        plt.plot(hds[:, int(z), 0, int(x)], label=name)
    plt.xlabel('Time')
    plt.ylabel('Heads [cm]')
    plt.legend(loc='lower right')

    # Plot concentrations
    plt.figure(2)
    for z, x, name in zip(z_bin, x_bin, coords.loc['well']['id']):
        plt.plot(ucn[:, int(z), 0, int(x)], label=name)
    plt.xlabel('Time')
    plt.ylabel('Concentrations [?]')
    plt.legend(loc='lower right')

    if show_plots:
        plt.show()

    print('Function execution finished.')

# Example usage:
# plot_heads_and_concentrations(scenario=2)