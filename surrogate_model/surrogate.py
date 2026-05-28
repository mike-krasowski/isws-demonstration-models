#%% IMPORTS ========================================
import flopy
import numpy as np
import matplotlib as mp; mp.use('TkAgg')
import matplotlib.pyplot as plt
import os
import pandas as pd

# import and run the toy model
from toy_model import run_toy_model
mf, run_time = run_toy_model()
# plt.close("all")

# quality of life
def vdir(var):
    for i in dir(var):
        print(i)

# access the head information from the model
import flopy.utils.binaryfile as bf
hds = bf.HeadFile(os.path.join(mf.model_ws, mf.name+'.hds')).get_alldata()

def xyzcellcenters_3d(mf):
    """
    Calculate 3D np.arrays of x-, y-, z- coodinate locations for each model cell node

    Parameters
    ----------
    mf - MODFLOW-NWT model object

    Returns
    -------
    xmesh - np.array of shape (layers, rows, columns), x-coordinate locations of model cell nodes

    ymesh - np.array of shape (layers, rows, columns), y-coordinate locations of model cell nodes

    zzz - np.array of shape (layers, rows, columns), z-elevations of model cell nodes
    """
    xxx, yyy, zzz = mf.modelgrid.xyzcellcenters
    xmesh = np.ones(zzz.shape)*-999
    ymesh = np.ones(zzz.shape)*-999
    for l in range(zzz.shape[0]):
        xmesh[l] = xxx
        ymesh[l] = yyy

    return xmesh, ymesh, zzz

def find_min_dist_and_value(mf, pkg_spd, key='stage'):
    """
    Calculates distance of each model cell to the nearest instance of a particular package cell (e.g., wel,
    riv) in the model. Also, returns the pertinent values (e.g., 'stage',  'flux') for that nearest package cell.

    Parameters
    ----------
    mf - MODFLOW-NWT model object

    pkg_spd - np.recarray, stress period data for a specified MODFLOW-NWT package (e.g., RIV, WEL, etc.)

    key - str, field/key for the np.recarray (`pkg_spd`) which denotes the pertinent simulated data to obtain from
    the model object.

    Returns
    -------
    dist - np.array of shape (layers, rows, columns), minimum distance from each cell to the nearest cell using the
    desired MF package (e.g., wel, riv)

    value - np.array of shape (layers, rows, columns), the pertinent (i.e., coordinated by key) value of the nearest
    package cell to each model cell

    """
    # package locations
    idx = (pkg_spd['k'], pkg_spd['i'], pkg_spd['j'])
    # model cell centers in 3d arrays
    xxx, yyy, zzz = xyzcellcenters_3d(mf)

    # find distance from package points
    stor1 = np.zeros((len(idx[0]), xxx.shape[0], xxx.shape[1], xxx.shape[2]))
    for n, i in enumerate(zip(idx[0], idx[1], idx[2])):
        stor1[n] = np.sqrt(((xxx - xxx[i]) ** 2) +
                           ((yyy - yyy[i]) ** 2) +
                           ((zzz - zzz[i]) ** 2))
    # calculate minimum distance and the values associated with package for each nearby cell
    dist = np.min(stor1, axis=0)
    min_idx = np.argmin(stor1, axis=0)
    value = pkg_spd[key][min_idx]
    return dist, value


# obtain 3D arrays of cell node xyz locations
xxx, yyy, zzz = xyzcellcenters_3d(mf)

for t in range(mf.nper):
    ### WELL DATA INTERPRETATIONS
    # well locations
    wel = mf.wel.stress_period_data[t]
    wel_dist, wel_Q = find_min_dist_and_value(mf, wel, key='flux')

    ### RIV CELL DATA INTERPRETATIONS
    # find riv cell locations
    riv = mf.riv.stress_period_data[t]
    riv_dist, riv_stage = find_min_dist_and_value(mf, riv, key='stage')

    ### RECHARGE INFORMATION
    rch = mf.rch.rech.array[t]
    rch_stor = np.zeros((xxx.shape))
    rch_stor[0,:] = rch[0]

    # create data frame of the current stress period
    data = pd.DataFrame(
                {'time': np.ones((xxx.size,))*t,
                'X': xxx.reshape(-1),
                'Y': yyy.reshape(-1),
                'Z': zzz.reshape(-1),
                # 'K_h': mf.upw.hk.array.reshape(-1),
                # 'K_v': mf.upw.vka.array.reshape(-1),
                # 'Ss': mf.upw.ss.array.reshape(-1),
                # 'Sy': mf.upw.sy.array.reshape(-1),
                # 'Rech': rch_stor.reshape(-1),
                'well_dist': wel_dist.reshape(-1),
                'well_flux': wel_Q.reshape(-1),
                'sw_dist': riv_dist.reshape(-1),
                'sw_stage': riv_stage.reshape(-1),
                'heads': hds[t,...].reshape(-1)}
        )
    if t == 0:
        df = data
    else:
        # concatenate dataFrame for additional stress periods
        df = pd.concat([df, data])

df.reset_index(drop=True, inplace=True)
mra = df.to_numpy() # MODFLOW results array

# keep only the heads of "active" cells
mra = mra[np.where(mra[:,-1] >= -999)]



# ======================================================================================================
### DIVIDE MODEL RESULTS RANDOMLY INTO TRAINING AND VALIDATION SETS
# create indices for splitting
rrr = np.random.rand(mra.shape[0], mra.shape[1])
split_idx = np.random.permutation(rrr.shape[0])
split_80  = int(mra.shape[0]*0.8)
# create training and test indices
train_idx, test_idx = split_idx[:split_80], split_idx[split_80:]
# split off training and test datasets
train_x, train_y = mra[train_idx,:-1], mra[train_idx, -1]
test_x,  test_y  = mra[test_idx,:-1],  mra[test_idx, -1]



# ======================================================================================================
### create ANN
import tensorflow as tf

ann = tf.keras.Sequential()

#AEJ is making this up - needs to read more on suggested model structure
ann.add(tf.keras.Input(shape=(train_x.shape[1],)))
for dl in [4,2,1,0.5]:
    ann.add(tf.keras.layers.Dense(int(train_x.shape[1]*dl))) # number of variables
ann.add(tf.keras.layers.Dense(1)) # matches number of output variables
ann.build()

ann.compile(optimizer='Adam', loss='mse')
print(
f'''
**ISWS Comment**
Model notes:
The model includes 5 layers that are shaped to mimic the number of output variables (i.e., 8) as defined
in the training dataset. Thus, the layer shapes are some multiple of 8 nodes - from 32 down to 1. 
(TBH, AEJ is not sure if this is correct, but this is where AEJ started.) There seems to be a link between
the number of parameters and the product of the node shapes of the previous and current layer: 

i.e., param_i = (layer_i.shape * layer_i-1.shape) + layer_i.shape 

However, again, more digging needs to be done to understand why this is the case.
**ISWS Comment**
'''
)
# print model summary
print(ann.summary())

# fit and evaluate the model
test = tf.keras.callbacks.EarlyStopping(
    monitor='loss',
    min_delta=1e-8,
    patience=0,
    verbose=0,
    mode='auto',
    baseline=None,
    restore_best_weights=False,
    start_from_epoch=0
)

import time
start = time.time()
ann.fit(train_x, train_y, epochs=30, callbacks=test)
print(f'Training time elapsed {start-time.time()}')
print('****************\n'+'ANN Evaluation...')
ann.evaluate(test_x, test_y)



################################################################################################################
## PLOTTING AND ANALYSIS

# plot the loss across the epochs
fig, ax = plt.subplots(1,1)
ax.plot(np.arange(len(ann.history.history['loss']))+1, ann.history.history['loss'])
ax.plot(np.arange(len(ann.history.history['loss']))+1, [0.001]*len(ann.history.history['loss']), '--r')
ax.set_ylabel('MSE')
ax.set_yscale('log')
ax.set_xlabel('epoch')
plt.grid('on')

# plot ann-predicted vs. modeled heads
ann_hds = ann.predict(mra[:,:-1])

for sp in np.unique(mra[:,0]):
    fig, ax = plt.subplots(3,3)
    for pl, z_depth in enumerate(np.unique(mra[:,3])[::-1]):
        subidx = np.where( (mra[:, 0] == sp) & (mra[:,3] == z_depth) )
        #plot model heads
        ccc = ax[pl,0].scatter(mra[subidx,1], mra[subidx,2], c=mra[subidx,-1])
        ax[pl,0].set_title(f'MF Heads:\nlayer {pl}, sp {sp}')
        cbar0 = fig.colorbar(ccc, ax=ax[pl,0])
        cbar0.set_label("MODFLOW simulated\nGW Head (ft)")

        # plot ann-heads
        ccc = ax[pl,1].scatter(mra[subidx, 1], mra[subidx, 2], c=ann_hds[subidx])
        ax[pl,1].set_title(f'ANN Heads:\nlayer {pl}, sp {sp}')
        ccc.set_clim([cbar0.vmin, cbar0.vmax])
        cbar1 = fig.colorbar(ccc, ax=ax[pl,1])
        cbar1.set_label("ANN simulated\nGW Head (ft)")
        # plot difference
        diff = mra[subidx,-1].reshape(-1)-ann_hds[subidx].reshape(-1)
        ccc = ax[pl,2].scatter(mra[subidx, 1], mra[subidx, 2],
                               c=diff,
                               cmap='bwr')
        ax[pl, 2].set_title(f'MF - ANN Heads:\nlayer {pl}, sp {sp}')
        ccc.set_clim([-0.1,0.1])
        cbar2 = fig.colorbar(ccc, ax=ax[pl, 2])
        cbar2.set_label("Heads Diff\nMF - ANN = Diff (ft)")
    # add space in subplots
    fig.subplots_adjust(hspace=0.45, wspace=0.25)




# Big Picture Next Steps:
# 1. get ANN to reproduce model results
# 2. can we modify inputs to ANN to predict model results?