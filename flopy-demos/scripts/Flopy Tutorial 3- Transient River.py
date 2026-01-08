#import packages
import flopy
import numpy as np
import matplotlib as mp; mp.use("TkAgg")
import matplotlib.pyplot as plt
import pandas as pd

#create model object
modelname = "my_river_model"

# path where files created for model run are to be saved
save_path  = ('./outputs/3_Transient_River/')
# points to the directory where the MODFLOW 2005 executable sits
# exe_path   = ('./modflowdir/mf2005.exe')
exe_path  = ("../../bin/win/mf2005.exe")

m = flopy.modflow.Modflow(modelname,
                          exe_name = exe_path,
                          model_ws = save_path )

#assign discretization variables
Lx = 1000.
Ly = 1000.
ztop = 0.
zbot = -50.
nlay = 1
nrow = 25
ncol = 25
dx = Lx/ncol
dy = Ly/nrow
dz = (ztop - zbot) / nlay

#specify number of stress periods
nper = 3

#specify if stress period is transient or steady-state
steady = [True, False, False] #<-- relates to [SS, transient transient]
perlen = [1.0, 5.0, 5.0]
nstp   = [1,   10,  10]

#create flopy discretization object, length units are meters (2) and time units are days (4)
dis = flopy.modflow.ModflowDis(model=m, nlay=nlay, nrow=nrow, ncol=ncol,
                               delr=dx, delc=dy, top=ztop, botm=zbot,
                               itmuni = 4, lenuni = 2,
                               nper=nper, steady=steady, perlen=perlen, nstp=nstp)

#CHECK GRID from DIS package
modelview = flopy.plot.PlotMapView(model=m, layer=0)
grid      = modelview.plot_grid()
# labeling the plot
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Grid', fontsize = 15, fontweight = 'bold')
plt.show()

#create ibound as array of ints = 1
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)

#create strt as array of floats = 1m
strt = np.ones((nlay, nrow, ncol), dtype=np.float32)

#create flopy bas object
bas = flopy.modflow.ModflowBas(m, ibound=ibound, strt=strt)

#define layer properties
laytyp = 1 #<-- a "convertible" layer
hk = 1.0 #m/d
vka = 1.0 #m/d
sy = 0.1 #1/m
ss = 1.e-4 #1/m

#attach LPF package
lpf = flopy.modflow.ModflowLpf(model=m, hk=hk, vka=vka, sy=sy, ss=ss, laytyp=laytyp, ipakcb=53)

#assign heads at start and end of stress period
strt_head=2
end_head=2
#create list to hold stress period constant head boundary condition cells
bound_sp1 = []

#assign constant head boundary cells on the left and right boundaries
for lay in range(nlay):
    for row in range(nrow):
        bound_sp1.append([lay,row,0,strt_head,end_head])
        bound_sp1.append([lay,row,ncol-1,strt_head,end_head])

print('Stress Period 1 CHD cells: \n',bound_sp1)

#create dictionary with stress period data
chd_spd={0: bound_sp1}

print('CHD Stress Period Data: \n', chd_spd)

#create flopy CHD object, and attach to model
chd = flopy.modflow.ModflowChd(model=m, stress_period_data=chd_spd)

#DEFINE RIVERS

#stress period 1 river cells
riv_sp1   = [] #create list to store all river cells for stress period 1
k_rivbott = 1 #river bottom hydraulic conductivity in m/d
sed_thick = 1 #thickness of riverbed sediment in m
cond      = k_rivbott*(dy)*(dx)/(sed_thick) #river bed conductance in m^2/d
r_stage   = 1 #stage in river (stress period 1)
r_bott    = 0 #river bottom
#assign river data to cells in central column
for iii in range(nrow):
    riv_sp1.append([0, iii, ncol/2, r_stage, cond, r_bott])


#stress period 2 river cells
riv_sp2 = [] #create list to store all river cells for stress period 2
r_stage = 5 #stage in river (stress period 2)
for iii in range(nrow):
    riv_sp2.append([0, iii, ncol/2, r_stage, cond, r_bott])

#stress period 3 river cells
riv_sp3 = [] #create list to store all river cells for stress period 3
r_stage = 1 #stage in river (stress period 3)
for iii in range(nrow):
    riv_sp3.append([0, iii, ncol/2, r_stage, cond, r_bott])


#create dictionary of stress period data
riv_spd = {0: riv_sp1, 1: riv_sp2, 2: riv_sp3}

#attach river package
riv = flopy.modflow.ModflowRiv(model=m,stress_period_data = riv_spd)

#CHECK BOUNDARY CONDITIONS
#use flopy to plot grid, ibound, rivers, and general head boundaries
modelview = flopy.plot.PlotMapView(model=m, layer=0)
grid      = modelview.plot_grid()
ib        = modelview.plot_ibound()
riv_plot  = modelview.plot_bc(ftype='RIV')
chd_plot  = modelview.plot_bc(ftype='CHD')
#add labels and legend
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Boundary Conditions', fontsize = 15, fontweight = 'bold')
plt.legend(handles=[mp.patches.Patch(color='green',label='River',ec='black'),
                    mp.patches.Patch(color='navy',label='Const Head Boundary',ec='black'),
                    mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                    mp.patches.Patch(color='black',label='Inactive Cell',ec='black')],
                    bbox_to_anchor=(1.8,1.0))
plt.show()

#create OC stress period data
oc_spd = {}
# creating dictionary of output instructions
for kper in range(nper):
    for kstp in range(nstp[kper]):
        oc_spd[(kper, kstp)] = ['save head',
                                'save drawdown',
                                'save budget',
                                'print head',
                                'print budget']
#creating output object
oc = flopy.modflow.ModflowOc(model=m,
                             stress_period_data=oc_spd,
                             compact=True)
# updating user on output instructions
print('Output Control \n', oc_spd)

#assign groundwater flow solver
pcg = flopy.modflow.ModflowPcg(model=m)

#write MODFLOW input files
m.write_input()

# Run the model
success, mfoutput = m.run_model(pause=False, report=True)
if not success:
    raise Exception('MODFLOW did not terminate normally.')

#extract binary data from head file
times   = [perlen[0],perlen[0]+perlen[1],perlen[0]+perlen[1]+perlen[2]] #extract times at end of each stress period
head    = {} #create dictionary to store head data at end of each stress period
frf     = {} #create dictionary to store flows through right cell face at end of each stress period
fff     = {} #create dictionary to store flows through front cell face at end of each stress period
headobj = flopy.utils.binaryfile.HeadFile(save_path+modelname+'.hds') #get head data as python object
budgobj = flopy.utils.binaryfile.CellBudgetFile(save_path+modelname+'.cbc') #get flow data as python object

#get data from python objects
for stress_per, time in enumerate(times): #iterate through times at end of each stress period
    head['sp%s'%(stress_per)] = headobj.get_data(totim=time) #append heads to head list for ea stress per
    frf['sp%s'%(stress_per)] = budgobj.get_data(text='FLOW RIGHT FACE',totim=time) #append right face flow to frf list for ea stress per
    fff['sp%s'%(stress_per)] = budgobj.get_data(text='FLOW FRONT FACE',totim=time) #append front face flow to fff list for ea stress per


# plot results for all stress periods
for i in range(len(times)):
    plt.figure(figsize=(9, 9))  # create 10 x 10 figure
    modelview = flopy.plot.PlotMapView(model=m, layer=0)  # use modelmap to attach plot to model
    grid = modelview.plot_grid()  # plot model grid
    riv_plot = modelview.plot_bc(ftype='RIV')  # plot river cells
    chd_plot = modelview.plot_bc(ftype='CHD')  # plot ghb cells
    contour_levels = np.linspace(head['sp%s' % i][0].min(), head['sp%s' % i][0].max(),
                                 15)  # set contour levels for contouring head
    head_contours = modelview.contour_array(head['sp%s' % i][0], levels=contour_levels)  # create head contours
    # flows          = modelview.plot_discharge(frf['sp%s'%i][0], fff['sp%s'%i][0], head=head['sp%s'%i]) #create discharge arrows

    # display parameters
    plt.xlabel('Lx (m)', fontsize=14)
    plt.ylabel('Ly (m)', fontsize=14)
    plt.title('Head/Flow Results: Stress Period %s' % (i + 1), fontsize=13, fontweight='bold')
    plt.colorbar(head_contours, aspect=5)
    plt.legend(handles=[mp.patches.Patch(color='green', label='River', ec='black'),
                        mp.patches.Patch(color='navy', label='Constant Head Boundary', ec='black'),
                        mp.patches.Patch(color='white', label='Active Cell', ec='black')],
               bbox_to_anchor=(1.8, 1.0))  # create legend
    plt.show()

#import 3d axes toolkit from matplotlib
from mpl_toolkits.mplot3d import Axes3D

#set X, Y, Z variables for 3d plot to be our model domain and head solution
X = np.arange(0,Lx,dx)
Y = np.arange(0,Ly,dy)
X, Y = np.meshgrid(X, Y)

#create a figure for every time
for i in range(len(times)):
    #create 3d figure
    fig_3d = plt.figure().add_subplot(projection='3d')
    Z = np.flipud(head['sp%s'%i][0]) #head matrix is flipped to display properly

    #create surface and labels
    surf = fig_3d.plot_surface(X,Y,Z, cmap = plt.cm.coolwarm, linewidth=0, antialiased=False, label='head') #plot head surface
    plt.colorbar(surf,shrink=0.5,aspect=5).set_label('Head (m)',fontsize=10,fontweight='bold') #set colorbar
    fig_3d.set_xlabel('Lx (m)', fontsize=15, fontweight='bold')
    fig_3d.set_ylabel('Ly (m)', fontsize=15, fontweight='bold')
    fig_3d.set_title('Head Surface: Stress-Period %s'%(i+1), fontsize=15, fontweight='bold')
    plt.show()

#plot a time series at cell left of river
#get time series for cell
cell_id = (0, int(nrow/2) - 1, int(ncol/2) - 1) #specify which cell we're interested in (cell at center of model)
time_series = headobj.get_ts(cell_id) #get the time series using flopy

#create plot
plt.figure()
plt.subplot(1, 1, 1)
plt.title('Head at cell ({0},{1},{2})'.format(cell_id[0] + 1,
                                              cell_id[1] + 1,
                                              cell_id[2] + 1),fontweight='bold') #create title with cell_id format
plt.xlabel('time (days)',fontweight='bold')
plt.ylabel('head (m)',fontweight='bold')
plt.plot(time_series[:, 0], time_series[:, 1], 'bo-') #plot the time series with points at each record
plt.show()
print('Time Series Head Data: \n', time_series) #print the time series data