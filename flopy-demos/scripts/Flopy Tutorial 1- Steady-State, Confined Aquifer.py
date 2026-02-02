import flopy
import numpy as np
import matplotlib as mp; mp.use("TkAgg")
import matplotlib.pyplot as plt


modelname = "my_model"
save_pth  = ('./outputs/1_SS_Confined/')
#On a different computer, update this directory pathway depending on the file location
# exe_path  = ("./modflowdir/mf2005.exe")
exe_path  = ("../../bin/win/mf2005.exe")

m = flopy.modflow.Modflow(modelname,
                          exe_name = exe_path,
                          model_ws = save_pth)

#assign discretization variables
Lx = 100.
Ly = 100.
ztop = 0.
zbot = -50.
nlay = 1     # number of layers
nrow = 10    # number of rows
ncol = 10    # number of columns
dx = Lx/ncol
dy = Ly/nrow
dz = (ztop - zbot) / nlay

#specify number of stress periods
nper = 1 # number of periods

#specify if stress period is transient or steady-state
steady = [True]
print("Steady-state data: \n", steady)

#create flopy discretization object, length and time are meters (2) and days (4)
dis = flopy.modflow.ModflowDis(model=m, nlay=nlay, nrow=nrow, ncol=ncol,
                               delr=dx, delc=dy, top=ztop, botm=zbot,
                               itmuni = 4, lenuni = 2,
                               nper=nper, steady=steady)

#use flopy to plot the grid of model 'm'
modelmap = flopy.plot.PlotMapView(model=m, layer=0)
grid = modelmap.plot_grid()
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Grid', fontsize = 15, fontweight = 'bold')
plt.show()

#create ibound as array of ints = 1
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
#assign left and right boundary cells to constant head
ibound[:, :, 0] = -1
ibound[:, :, -1] = -1

print("ibound values: \n", ibound)

#create strt as array of floats = 1
strt = 5*np.ones((nlay, nrow, ncol), dtype=np.float32)
#set left side head to 10 m
strt[:, :, 0] = 10.
#set right side head to 0 m
strt[:, :, -1] = 0.

print("starting head values: \n", strt)

#create flopy bas object
bas = flopy.modflow.ModflowBas(m, ibound=ibound, strt=strt)

#plot grid and ibound
modelmap = flopy.plot.PlotMapView(model=m, layer=0)
grid = modelmap.plot_grid()
ib = modelmap.plot_ibound()
#add labels and legend
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Ibound', fontsize = 15, fontweight = 'bold')
plt.legend(handles=[mp.patches.Patch(color='blue',label='Const. Head',ec='black'),
                   mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                   mp.patches.Patch(color='black',label='Inactive Cell',ec='black')],
                   bbox_to_anchor=(1.5,1.0))
plt.show()

#define horizontal hydraulic conductivity
hk = np.ones((nlay,nrow,ncol), dtype=np.float32)
vk = np.ones((nlay,nrow,ncol), dtype=np.float32)
print("horizontal k values: \n", hk,
     "\n vertical k values: \n", vk)

#define layer type as confined
laytyp = np.zeros((nlay,), dtype=np.int32)
print("layer type values: \n", laytyp)

lpf = flopy.modflow.ModflowLpf(model=m, hk=hk, vka=vk, laytyp=laytyp, ipakcb=1)

#create oc stress period data.
spd = {(0, 0): ['save head', 'save budget']}
print("oc stress period data: \n", spd)

oc = flopy.modflow.ModflowOc(model=m, stress_period_data=spd, compact=True)

pcg = flopy.modflow.ModflowPcg(model=m)

#write MODFLOW input files
m.write_input()
# Run the model
success, mfoutput = m.run_model(pause=False, report=True)
if not success:
    raise Exception('MODFLOW did not terminate normally.')

#extract binary data from head file as flopy head object
headobj = flopy.utils.binaryfile.HeadFile(save_pth + modelname+'.hds')
print("flopy head object: \n", headobj)

#extract head data from head object
head = headobj.get_data(totim=1.0)

#extract binary data from budget file as flopy budget object
budgobj = flopy.utils.binaryfile.CellBudgetFile(save_pth+modelname+'.cbc')
print("flopy budget object: \n", budgobj)

frf = budgobj.get_data(text='flow right face', totim=1.0)
fff = budgobj.get_data(text='flow front face', totim=1.0)

print("Flow through Right Face of Grid Cells m^3/d \n", frf,
     "\n Flow through Front Face of Grid Cells m^3/d \n", fff)

#plot results
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelmap = flopy.plot.PlotMapView(model=m, layer=0) #use modelmap to attach plot to model
grid = modelmap.plot_grid() #plot model grid
contour_levels = np.linspace(head[0].min(),head[0].max(),11) #set contour levels for contouring head
head_contours = modelmap.contour_array(head, levels=contour_levels) #create head contours

#display parameters
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Steady-State Model, Flow(m^3/d) and Head(m) Results', fontsize = 15, fontweight = 'bold')
plt.colorbar(head_contours,aspect=5)
plt.show()

#import 3d axes toolkit from matplotlib
from mpl_toolkits.mplot3d import Axes3D

#create 3d figure
fig_3d = plt.figure().add_subplot(projection='3d')

#set X, Y, Z variables for 3d plot to be our model domain and head solution
X = np.arange(0,Lx,dx)
Y = np.arange(0,Ly,dy)
X, Y = np.meshgrid(X, Y)
Z = np.flipud(head[0])

#create surface and labels
surf = fig_3d.plot_surface(X,Y,Z, cmap = 'viridis', linewidth=0, label='head')
plt.colorbar(surf,shrink=0.5,aspect=5).set_label('Head (m)',fontsize=10,fontweight='bold')
fig_3d.set_xlabel('Lx (m)', fontsize=15, fontweight='bold')
fig_3d.set_ylabel('Ly (m)', fontsize=15, fontweight='bold')
fig_3d.set_title('Steady-State Model Head Profile', fontsize=15, fontweight='bold')
plt.show()