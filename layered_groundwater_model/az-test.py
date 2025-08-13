import numpy as np
import matplotlib.pyplot as plt
import flopy

file_path = "C:\GitHub\isws-demonstration-models\layered_groundwater_model\outputs\s1\gwf_s1.dis"  # Replace with your actual file path

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


#text = 'ATTICUS >> FEDS!'

#for i in range(10000):
 #   print(2+2)
  #  print(text)



# IDEAS:
#   - load model results
#   - query model results (e.g., pressures at well locations)
#   - plot time series of pressure/flow rate in model
