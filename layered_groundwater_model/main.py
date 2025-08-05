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

    start_time = time.time()

    success = False
    wiggle = 0
    while not success:

        success, sim, nodes = run_the_models(name, wiggle)
        wiggle = random.choice([1,-1]) * random.random() * 0.1

    sim = make_the_animation(sim, nodes)

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

# user inputs, user inputs, user inputs
max_workers = 9

spds_path = './inputs/scenarios_short3.xlsx'
exfi = pd.ExcelFile(spds_path)
scenarios = exfi.sheet_names

scenarios = scenarios[:2]

# end user inputs, end user inputs, end user inputs

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


