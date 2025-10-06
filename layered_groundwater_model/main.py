"""
The purpose of this script is to build MODFLOW models to simulate the various scenarios of the table-top layered
groundwater model and produce animations of the results

author: ISWS, krasows2, zavelle
July 14th, 2025
"""

# imports
import concurrent.futures
import flopy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
import pandas as pd
import random
from shapely import Polygon
import shutil
import time

from modules import *

def main(name):
    """
    The purpose of this function is to coordinate running an individual scenario and constructing its animation in a
    single function that can be parallelized. The "wiggle" is being applied because there is some convergence
    sensitivity to small changes in input parameters. "Wiggling" model parameters ensures we get a solution for each
    scenario. This does mean the parameters between scenarios aren't necessarily identical.
    :param name: str
        the name of the scenario. this will also be use as the simulation/model name
    :return: list
        two element list containing the simulation in position 0 and the runtime in position 1
    """

    start_time = time.time()

    success = False
    wiggle = 0
    while not success:

        success, sim, nodes = run_the_models(name, wiggle)
        wiggle = random.choice([1,-1]) * random.random() * 0.1

    sim = make_the_animation(sim, nodes, wiggle, parameter='concentration')

    return [sim, time.time() - start_time]

# housekeeping
start_time = time.time()

# define an outputs folder
savepath = 'outputs'

# if it exists, clear it
if os.path.exists(savepath):
    shutil.rmtree(savepath)

# make the outputs folder
os.makedirs(savepath)

## =================================================================================
## COMMON USER INPUTS - COMMON USER INPUTS - COMMON USER INPUTS - COMMON USER INPUTS
## =================================================================================

# how many cores to work in parallel
max_workers = 9

# which spreadsheet we would like to read our scenarios from
spds_path = 'inputs/scenarios.xlsx'

## =================================================================================
## END COMMON USER INPUTS
## =================================================================================

# create excel file object to get the sheet names to then look at each sheet
exfi = pd.ExcelFile(spds_path)
scenarios = exfi.sheet_names

# only some scenarios are constructed for now so only run those
scenarios = scenarios[:3]

if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(main, scenarios)

func_time = 0
sims = []
for ppp, result in enumerate(results):

    sims.append(result[0])
    func_time += result[1]

print('ISWS: the WHOLE script took: {} mins with {} mins of function time'.format(
    round((time.time() - start_time)/60, 2),
    round(func_time/60, 2)
))


