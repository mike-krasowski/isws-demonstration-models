
"""The purpose of this script is to implement contaminant transport (saturated zone) into the ParFlow model.

It is created to be very similar to the ISWS State Fair model. The hydrological units in the ParFlow contaminant transport model are the same as those in the ISWS State Fair model. 
This model is built using ParFlow version 3.13.0.To run this contamination model, please first create the PFB input folder using create_pfb.ipynb script.

Contamination_model.ipynb contains scripts that create ParFlow models and animations representing the enVISION tabletop groundwater
model. The workflow for running the ParFlow model(s) and producing the animations currently relies on using Docker due 
to ParFlow's requirement for running on linux.

September 15, 2025

Author(s): ISWS; alljones, esrag. """