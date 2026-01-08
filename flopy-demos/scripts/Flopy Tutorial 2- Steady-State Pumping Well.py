#import packages
import flopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp

#create model object
modelname  = "my_pumping_model"

# path where files created for model run are to be saved
save_path  = ('./outputs/2_Pumping_SS/')
# points to the directory where the MODFLOW 2005 executable sits
# exe_path   = ('./modflowdir/mf2005.exe')
exe_path  = ("../../bin/win/mf2005.exe")

m = flopy.modflow.Modflow(modelname,
                          exe_name = exe_path,
                          model_ws = save_path)

#assign discretization variables
Lx = 100.
Ly = 100.
ztop = 0.
zbot = -50.
nlay = 1
nrow = 25
ncol = 25
dx = Lx/ncol
dy = Ly/nrow
dz = (ztop - zbot) / nlay

#specify number of stress periods
nper = 1

#specify if stress period is transient (False) or steady-state (True)
steady = [True]

#create flopy discretization object, length and time are meters (2) and days (4)
dis = flopy.modflow.ModflowDis(model=m, nlay=nlay, nrow=nrow, ncol=ncol,
                               delr=dx, delc=dy, top=ztop, botm=zbot,
                               itmuni = 4, lenuni = 2,
                               nper=nper, steady=steady)

# Use PlotMapView instead
mapview = flopy.plot.map.PlotMapView(model=m)
grid = mapview.plot_grid()
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Grid', fontsize = 15, fontweight = 'bold')
plt.show()

#create ibound as array of ints = 1
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
#assign left and right boundary cells to constant head
ibound[:, :, 0] = -1
ibound[:, :, -1] = -1
# assign northern- and southern-most boundary cells to constant head
ibound[:,0,:] = -1
ibound[:,-1,:] = -1

#create strt as array of floats = 10m <-- initial conditions
strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
strt[:, :, :] = 10.

#create flopy bas object
bas = flopy.modflow.ModflowBas(m, ibound=ibound, strt=strt)

#CHECK IBOUND
#use flopy to plot grid and ibound
modelview = flopy.plot.PlotMapView(model=m, layer=0)
grid      = modelview.plot_grid()
ib        = modelview.plot_ibound()
#add labels and legend
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Ibound', fontsize = 15, fontweight = 'bold')
plt.legend(handles=[mp.patches.Patch(color='blue',label='Const Head',ec='black'),
                    mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                    mp.patches.Patch(color='black',label='Inactive Cell',ec='black')],
                    bbox_to_anchor=(1.5,1.0))
plt.show(modelview)

#define horizontal and vertical hydraulic conductivity
hk = np.ones((nlay,nrow,ncol), dtype=np.float32) # assumes horizontal isotropy
vk = np.ones((nlay,nrow,ncol), dtype=np.float32)

#define specific storage
ss = np.ones((nlay,nrow,ncol), dtype=np.float)
ss[:,:,:] = 1e-5

#define layer type as confined
laytyp = np.zeros((nlay,), dtype=np.int32)

#create flopy layer property flow object
lpf = flopy.modflow.ModflowLpf(model=m, hk=hk, vka=vk, ss=ss, laytyp=laytyp, ipakcb=1)

#Create Single Well at center of domain with [lay, row, col, flux] list
pumping_rate = -100 #in m^3/d, negative for pumping/positive for injection
well_1 = [0, ncol/2, nrow/2, pumping_rate]
print("Well 1 [layer, row, column, flux]: \n", well_1)

#Create Dictionary With Stress Period Data
wel_spd = {0: [well_1]}
print("Well Stress Period Data: \n", wel_spd)

#Create flopy wel object
wel = flopy.modflow.ModflowWel(model=m, stress_period_data=wel_spd)

#CHECK WELL LOCATION
#use flopy to plot grid, ibound, and wells
modelview = flopy.plot.PlotMapView(model=m, layer=0)
grid      = modelview.plot_grid()
ib        = modelview.plot_ibound()
wel       = modelview.plot_bc(ftype='WEL')

#add labels and legend
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Well and Boundary Conditions', fontsize = 15, fontweight = 'bold')
plt.legend(handles=[mp.patches.Patch(color='red',label='Well',ec='black'),
                    mp.patches.Patch(color='blue',label='Const Head',ec='black'),
                    mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                    mp.patches.Patch(color='black',label='Inactive Cell',ec='black')],
                    bbox_to_anchor=(1.5,1.0))
plt.show(modelview)

#create oc stress period data.
spd = {(0, 0): ['print head', 'print budget', 'save head', 'save budget']}
# above, we are specifying that for the (1st period and 1st timestep) the heads and budget are saved and printed

#create flopy output control object
oc = flopy.modflow.ModflowOc(model=m, stress_period_data=spd, compact=True)

#assign groundwater flow solver
pcg = flopy.modflow.ModflowPcg(model=m)

#write MODFLOW input files
m.write_input()

# Run the model
success, mfoutput = m.run_model(pause=False, report=True)
if not success:
    raise Exception('MODFLOW did not terminate normally.')

#extract binary data from head file
headobj = flopy.utils.binaryfile.HeadFile(save_path+modelname+'.hds')
head = headobj.get_data(totim=1.0)

#extract binary data from budget file
budgobj = flopy.utils.binaryfile.CellBudgetFile(save_path+modelname+'.cbc')
frf = budgobj.get_data(text='flow right face', totim=1.0)
fff = budgobj.get_data(text='flow front face', totim=1.0)

# initiating Latex rendering
# mp.rcParams['text.usetex'] = True

#plot results
plt.figure(figsize=(10,10)) #create 10 x 10 figure
modelview      = flopy.plot.PlotMapView(model=m, layer=0) #use modelmap to attach plot to model
grid           = modelview.plot_grid() #plot model grid
contour_levels = np.linspace(head[0].min(),head[0].max(),15) #set contour levels for contouring head
head_contours  = modelview.contour_array(head, levels=contour_levels) #create head contours
## flows          = modelview.plot_discharge(frf[0], fff[0], head=head) #create discharge arrows

#display parameters
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Steady-State Pumping, Flow(m^3/d) and Head(m) Results', fontsize = 15, fontweight = 'bold')
plt.colorbar(head_contours, aspect=5)
plt.show(modelview)

#create plot
fig, ax = plt.subplots(figsize=(9,7))
cf = plt.contourf(np.flipud(head[0,:,:]), levels=contour_levels)
#display parameters
plt.xlabel('Lx (m)',fontsize = 14) # <-- originally displays the column number
ax.set_xticklabels( ax.get_xticks()*dx )
plt.ylabel('Ly (m)',fontsize = 14) # <-- originally displays the row number
ax.set_yticklabels( ax.get_yticks()*dy )
plt.title('Steady-State Pumping, Head(m) Results', fontsize = 15, fontweight = 'bold')
plt.colorbar(cf,aspect=10)
plt.show()

#import 3d axes toolkit from matplotlib
from mpl_toolkits.mplot3d import Axes3D

#create 3d figure
fig_3d = plt.figure(figsize=(12,5))
ax = fig_3d.gca(projection='3d')

#set X, Y, Z variables for 3d plot to be our model domain and head solution
X = np.arange(0,Lx,dx)
Y = np.arange(0,Ly,dy)
X, Y = np.meshgrid(X, Y)
Z = np.flipud(head[0])

#create surface and labels
surf = ax.plot_surface(X,Y,Z, cmap = plt.cm.coolwarm, linewidth=0, antialiased=False, label='head')
fig_3d.colorbar(surf,shrink=0.5,aspect=5).set_label('Head (m)',fontsize=10,fontweight='bold')
ax.set_xlabel('Lx (m)', fontsize=15, fontweight='bold')
ax.set_ylabel('Ly (m)', fontsize=15, fontweight='bold')
ax.set_title('Steady-State Pumping, Head Surface', fontsize=15, fontweight='bold')
plt.show(surf)

#plot head head transect
plt.figure(figsize = (10,4))
x = np.arange(0,Lx,dx)
plt.plot(x,np.flipud(head[0])[int(nrow/2)][:])
plt.title('Head Transect across X-Domain',fontweight = 'bold', fontsize = 14)
plt.xlabel('X Distance (m)',fontsize = 12)
plt.ylabel('Head (m)',fontsize = 12)
plt.show()