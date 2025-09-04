#import numpy as np
import matplotlib as mpl; mpl.use("WebAgg")
import matplotlib.pyplot as plt
import flopy

def vdir(n):
    for i in dir(n):
        print(i)
#############pick scenario to run
scenario = 2

file_path = f"./outputs/s{scenario}/"  # Replace with your actual file path

mf = flopy.mf6.MFSimulation.load(
    sim_name="mfsim", #.nam",
    sim_ws=file_path,
    )
gwf = mf.get_model(f"gwf_s{scenario}")
gwt = mf.get_model(f"gwt_s{scenario}")

####just a note for atticus on syntax for pointing at files
# pathhds = r".\outputs\s2\gwf_s2.hds"


#import flopy.utils.binaryfile as bf
import pandas as pd
# hds = bf.HeadFile(pathhds).get_alldata() # (time, layer, row, column)
hds = gwf.output.head().get_alldata()
ucn = gwt.output.concentration().get_alldata()

print("grabbed data heads and concentration data ")

# read in well csv --> pandas as pd --> pd.read_csv
well = pd.read_csv('./inputs/well_info.csv')

# determine well locations and cells from csv
coords = well.set_index('type')
x_bin = coords.loc['well']['x_coord'] // gwf.modelgrid.delr[0]
z_bin = coords.loc['well']['z_coord'] // gwf.modelgrid.delz[0,0,0]

#####plot for Headsfile
# plot time series at each location --> plt.plot(hds[time, z, y(into page), x])
plt.figure(1)
for z, x, name in zip(z_bin, x_bin, coords.loc['well']['id']):
    plt.plot(hds[:, int(z), 0, int(x)], label=name)
plt.xlabel('Time')
plt.ylabel('Heads [cm]')
plt.legend(loc='lower right')

##########plot for UcnFIle
# plot time series at each location --> plt.plot(ucn[time, z, y(into page), x])
plt.figure(2)
for z, x, name in zip(z_bin, x_bin, coords.loc['well']['id']):
    plt.plot(ucn[:, int(z), 0, int(x)], label=name)
plt.xlabel('Time')
plt.ylabel('Concentrations [?]')
plt.legend(loc='lower right')


plt.show()
plt.close('all')
print('Script Finished.')


# except FileNotFoundError:
#     print(f"Error: The file at '{file_path}' was not found.")
# except Exception as e:
#     print(f"An error occurred: {e}")
