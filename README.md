# isws-demonstration-models
This repository is for computer models meant to replicate the behavior of physical demonstration models used by the 
Illinois State Water Survey. Each modeling effort will be housed in its own folder. As a result, each folder within this
repository will have its own _config_ folder with an _environment.yml_ file for setting the appropriate Python 
environment. See the subsections below for descriptions and instructions specific to each folder of this repo.

### bin
This is a shared folder for the storing of important binary files (such as MODFLOW executables) needed for the various 
folder of this repository.

### flopy-demos
The scripts and files of this folder contain examples that demonstrate the usage of FloPy and methods for investigating 
model outputs. 

Note, many of the scripts in this folder were assembled by Daniel Abrams and Shelby Ahrendt in 2019 and earlier and 
hence may require older versions of MODFLOW and FloPy to function. Future work may see these updated.

### layered_groundwater_model_modflow
This folder contains scripts that create MODFLOW 6 models and animations representing the enVISION tabletop groundwater
model. The model and scripts have infrastructure to build and a run groundwater flow, tranpsort, and energy models.

### layered_groundwater_model_parflow
This folder contains scripts that create ParFlow models and animations representing the enVISION tabletop groundwater
model. The workflow for running the ParFlow model(s) and producing the animations currently relies on using Docker due 
to ParFlow's requirement for running on linux.

