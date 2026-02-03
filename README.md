# isws-demonstration-models
This repository is for computer models meant to demonstrate the use of various modeling software used by the 
Illinois State Water Survey. Each modeling effort will be housed in its own folder. As a result, each folder within this
repository will have its own _config_ folder with an _environment.yml_ file for setting the appropriate Python 
environment. See the subsections below for descriptions and instructions specific to each folder of this repo.

### bin
This is a shared folder for the storing of important binary files (such as MODFLOW executables) needed for the various 
folder of this repository.

### layered_groundwater_model
This folder contains scripts that create MODFLOW 6 models and animations representing the enVISION tabletop groundwater
model. The model and scripts have infrastructure to build and a run groundwater flow, tranpsort, and energy models.
Please see the readme in this folder for more information.
