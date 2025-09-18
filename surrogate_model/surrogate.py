#%% IMPORTS ========================================
import flopy
import numpy as np
import matplotlib as mp; mp.use('TkAgg')
import matplotlib.pyplot as plt
import os

# import and run the toy model
from toy_model import run_toy_model
mf, run_time = run_toy_model()

# quality of life
def vdir( var ):
    for i in dir(var):
        print(i)

# access the head information from the model
import flopy.utils.binaryfile as bf
hds = bf.HeadFile(os.path.join(mf.model_ws, mf.name+'.hds')).get_alldata()

# Parameters to train the surrogate model on:
# >> time
# >> location - x, y, z
# >> distance to nearest well? --> nearest 5 wells?
# >> pumping rate of nearest well? --> nearest 5 wells?
# >> distance to nearest surface water?
# >> elevation of surface water?

# create dataframe of values
def xyzcellcenters_3d(mf):
    xxx, yyy, zzz = mf.modelgrid.xyzcellcenters
    xmesh = np.ones(zzz.shape)*-999
    ymesh = np.ones(zzz.shape)*-999
    for l in range(zzz.shape[0]):
        xmesh[l] = xxx
        ymesh[l] = yyy

    return xmesh, ymesh, zzz


def find_min_dist_and_value(mf, pkg_spd, key='stage'):

    idx = (pkg_spd['k'], pkg_spd['i'], pkg_spd['j'])
    xxx, yyy, zzz = xyzcellcenters_3d(mf)

    # find distance from surface waters
    stor1 = np.zeros((len(idx[0]), xxx.shape[0], xxx.shape[1], xxx.shape[2]))
    for n, i in enumerate(zip(idx[0], idx[1], idx[2])):
        stor1[n] = np.sqrt(((xxx - xxx[i]) ** 2) +
                           ((yyy - yyy[i]) ** 2) +
                           ((zzz - zzz[i]) ** 2))

    dist = np.min(stor1, axis=0)
    min_idx = np.argmin(stor1, axis=0)
    value = riv[key][min_idx]
    return dist, value



xxx, yyy, zzz = xyzcellcenters_3d(mf)

### WELL DATA INTERPRETATIONS
# well locations
wel = mf.wel.stress_period_data[0]
idx = (wel['k'],wel['i'],wel['j'])
wel_Q = wel['flux']

# all cells distance from well location
wel_dist = np.sqrt( ((xxx - xxx[idx])**2)+
                    ((yyy - yyy[idx])**2)+
                    ((zzz - zzz[idx])**2) )


### RIV CELL DATA INTERPRETATIONS
# find riv cell locations
riv = mf.riv.stress_period_data[0]
idx = (riv['k'], riv['i'], riv['j'])

# find distance from surface waters
stor1 = np.zeros((len(idx[0]), xxx.shape[0], xxx.shape[1], xxx.shape[2]))
for n, i in enumerate(zip(idx[0], idx[1], idx[2])):
    stor1[n] = np.sqrt( ((xxx - xxx[i])**2)+
                        ((yyy - yyy[i])**2)+
                        ((zzz - zzz[i])**2) )

riv_dista  = np.min(stor1, axis=0)
min_idx   = np.argmin(stor1, axis=0)
riv_stagea = riv['stage'][min_idx]

riv_distb, riv_stageb = find_min_dist_and_value(mf, riv, key='stage')

print( 'Dist: ', (riv_dista==riv_distb).all() )
print( 'Stage: ', (riv_stagea==riv_stageb).all() )


# Next steps:
# 1. function - nearest calculations (all wells/all rivs)
# 2. loop stress periods











