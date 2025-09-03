import numpy as np
import matplotlib as mpl; mpl.use("WebAgg")
import matplotlib.pyplot as plt
import flopy

def vdir(n):
    for i in dir(n):
        print(i)

file_path = r"./outputs/s1/"  # Replace with your actual file path

mf = flopy.mf6.MFSimulation.load(
    sim_name="mfsim", #.nam",
    sim_ws=file_path,
    )
gwf = mf.get_model("gwf_s1")

path2 = r".\outputs\s1\gwf_s1.hds"

import flopy.utils.binaryfile as bf
import pandas as pd
hds = bf.HeadFile(path2).get_alldata() # (time, layer, row, column)

print ("squeeks")

# read in well csv --> pandas as pd --> pd.read_csv
well = pd.read_csv('./inputs/well_info.csv')
print (well)



# # determine well locations and cells from csv
# xyz = gwf.modelgrid.xyzcellcenters
# coords = well.loc[
#     ((well['type']=='well') & (well['id']=='well_2')),
#     ['id','x_coord', 'z_coord']
# ]
# print(coords)
#
# x_bin = int(coords.iloc[0]['x_coord'] // gwf.modelgrid.delr[0])
# z_bin = int(coords.iloc[0]['z_coord'] // gwf.modelgrid.delz[0,0,0])
# print(x_bin)
# print(z_bin)
#
# # plot time series at each location --> plt.plot(hds[time, z, y(into page), x])
# plt.plot(hds[:, z_bin, 0, x_bin])
# plt.show()
#
# print('Script Finished.')

# determine well locations and cells from csv
xyz = gwf.modelgrid.xyzcellcenters
coords = well.set_index('type')

x_bin = coords.loc['well']['x_coord'] // gwf.modelgrid.delr[0]
z_bin = coords.loc['well']['z_coord'] // gwf.modelgrid.delz[0,0,0]
print(x_bin)
print(z_bin)

# plot time series at each location --> plt.plot(hds[time, z, y(into page), x])
plt.figure()
for z, x, name in zip(z_bin, x_bin, coords.loc['well']['id']):
    plt.plot(hds[:, int(z), 0, int(x)], label=name)
plt.legend(loc='lower right')
plt.show()

print('Script Finished.')




'''
try:
    with open(file_path, 'r') as file:
        # Read the entire content of the file
        content = file.read()
        print("File content:")
        print(content)

        # Or read line by line
        #for line in file:
         # print(line.strip()) # .strip() removes leading/trailing whitespace including newlines

except FileNotFoundError:
    print(f"Error: The file at '{file_path}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")
'''

#text = 'ATTICUS >> FEDS!'

#for i in range(10000):
 #   print(2+2)
  #  print(text)



# IDEAS:
#   - load model results
#   - query model results (e.g., pressures at well locations) (gwf)
#   - plot time series of pressure/flow rate in model (gwf)
#   - plot breakthrough curves (gwt)
