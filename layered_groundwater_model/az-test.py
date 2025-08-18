import numpy as np
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
print (hds)
print("soy")
df = pd.DataFrame(hds)
df.to_csv("heads_data.csv", index=None)


print ("squeeks")
# read in well csv --> pandas as pd --> pd.read_csv
# determine well locations and cells from csv
# plot time series at each location

# plt.plot(hds[time, z, y(into page), x])
#plt.plot(hds[:,30, 0, 50])
#plt.show()

raise Exception('AEJ breaks stuff.')

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
