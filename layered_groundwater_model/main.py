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

# housekeeping
start_time = time.time()

def main(input):
    """
    The purpose of this function is to coordinate running an individual scenario and constructing its animation in a
    single function that can be parallelized. The "wiggle" is being applied because there is some convergence
    sensitivity to small changes in input parameters. "Wiggling" model parameters ensures we get a solution for each
    scenario we submit. This means the parameters between scenarios aren't necessarily identical, but we're making a
    nice visualization, not a model of a real world problem.
    :param input: tuple
        two entry tuple with the first being a str indicating the scenario/model name and the second being a list of the
        output parameters to be animated (i.e. 'head', 'concentration', or 'temperature')
    :return: list
        two element list containing the simulation in position 0 and the runtime in position 1
    """

    # unpack the input tuple
    name = input[0]
    parameters = input[1]

    # for keeping track of how much time we are saving by parallelizing
    start_time = time.time()

    # run the model until we have a successful run
    success = False
    wiggle = 0
    while not success:

        success, sim, nodes = build_and_run_models(name, wiggle)
        wiggle = random.choice([1,-1]) * random.random() * 0.1

    # create animations for each parameter requested
    for parameter in parameters:
        make_animation(sim, nodes, wiggle, parameter=parameter)

    return [sim, time.time() - start_time]

## =================================================================================
## COMMON USER INPUTS - COMMON USER INPUTS - COMMON USER INPUTS - COMMON USER INPUTS
## =================================================================================

# how many cores to work in parallel
max_workers = 6

# which spreadsheets to read our scenarios from. The script will look for these csvs at ./inputs/scenarios
# 's1' - long-running flow demonstration,
# 's2' - injection at well 2, extraction at well 6
scenarios = ['s2']

# options here are 'head' (gwf), 'concentration' (gwt), and 'temperature' (gwe)
parameters = ['concentration']

## =================================================================================
## END COMMON USER INPUTS
## =================================================================================

# define an outputs folder
savepath = 'outputs'

# if it exists, remove it and all of its contents
if os.path.exists(savepath):
    shutil.rmtree(savepath)

# make the outputs folder
os.makedirs(savepath)

# if there's only one scenario, do not try to run in parallel
if len(scenarios) == 1:

    result = main((scenarios[0], parameters))
    results = [result]

else:

    # for now, we know there are only two scenarios established, run them both.
    inputs = [(scenarios[0], [None, 'temperature']),
              (scenarios[1], ['concentration', 'temperature'])]

    if __name__ == "__main__":
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = executor.map(main, inputs)

func_time = 0
sims = []
for ppp, result in enumerate(results):

    sims.append(result[0])
    func_time += result[1]

print('ISWS: the WHOLE script took: {} mins with {} mins of main function time'.format(
    round((time.time() - start_time)/60, 2),
    round(func_time/60, 2)
))
