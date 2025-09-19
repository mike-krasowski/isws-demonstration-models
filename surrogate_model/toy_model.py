"""
Simplified model to test the linkages between MT3D and a transient model.

Author: Allan Jones
"""

#%% IMPORTS ========================================
import flopy
import numpy as np
import matplotlib as mp; mp.use('TkAgg')
import matplotlib.pyplot as plt
import os

def vdir( var ):
    for i in dir(var):
        print(i)


#%% INITIALIZE IMPORTANT DIRECTORIES AND THE MODEL OBJECT ========================================
def run_toy_model(model_name=None):
    # define the model name
    if model_name == None:
        model_name = "toy-Model"
    # directory of MODFLOW executable
    exe_dir = '../bin/win/mfnwt.exe'

    # directory for saving model files
    from datetime import datetime
    folder_date   = datetime.now().strftime('%Y%m%d')
    modelfile_dir = (f'./outputs/{model_name}')

    # create the model object
    m = flopy.modflow.Modflow( model_name,
                               exe_name=exe_dir,
                               model_ws=modelfile_dir,
                               version="mfnwt",
                               stop_tol=1e-8 )
    # set units
    m.set_model_units(iunit0=1)

    #%% ASIGN NWT AS SOLVER ========================================
    #assign groundwater flow solver
    nwt = flopy.modflow.ModflowNwt( m,
                                    iprnwt=1,
                                    maxiterout=1000,
                                    linmeth=2,
                                    )

    #%% INITIALIZE THE GRID AND DISCRETIZE MODEL - SIMPLE ========================================
    #assign discretization variables
    Lx = 100. #x dimension (m)
    Ly = 100. #y dimension (m)
    ztop = 0. #top of model (m)
    zbot = -12. #bottom of model (m)
    nlay = 3 #number of model layers
    nrow = 25 #number model rows (cuts up model along Ly dimension)
    ncol = 25 #number of model columns (cuts up model along Lx dimension)
    dx = Lx/ncol #grid spacing along Lx dimension
    dy = Ly/nrow #grid spacing along Ly dimension
    dz = (ztop - zbot) / nlay #grid spacing between layers

    #specify number of stress periods
    nper = 5 #!!! eventually update to include transient

    #specify if stress period is transient (False) or steady-state (True)
    steady = []
    for i in range(nper):
        if i > 0:
            steady.append(False)
        else:
            steady.append(True)

    nstp = list(np.ones(nper,  dtype=int) * 365)
    nstp[0] = 1

    #create flopy discretization object, length and time are meters (2) and days (4)
    model_btm = np.ones((nlay,nrow,ncol))
    for i in range(nlay):
        model_btm[i] = model_btm[i]*(ztop-((i+1)*dz))
    dis = flopy.modflow.ModflowDis(model=m,
                                   nlay=nlay,
                                   nrow=nrow,
                                   ncol=ncol,
                                   delr=dx,
                                   delc=dy,
                                   top=ztop,
                                   botm=model_btm,
                                   itmuni = 4,
                                   lenuni = 2,
                                   nper=nper,
                                   steady=steady)

    #%% CREATE OUTPUT CONTROL (OC) PACKAGE  ========================================
    oc_spd = {} #create oc stress period data.
    for kper in range(nper):
        for kstp in range(nstp[kper]):
            oc_spd[(kper, kstp)] = ['save head',
                                    'save drawdown',
                                    'save budget',
                                    'print head',
                                    'print budget' ]
    #create flopy output control object
    oc = flopy.modflow.ModflowOc(model=m, stress_period_data=oc_spd, compact=True)


    #%% ESTABLISH RIVER CELLS - WESTERN BOUNDARY ========================================
    #create list to store river stress period data
    stage_riv = np.cos( np.array(range(nper))*np.pi/2 ) + 7 # set west river stage to oscillate around 5 m
    lake_stage = np.ones((nper,))*5.5 #np.sin( np.array(range(nper))*np.pi/2 ) + 5.5
    cond = 1  # set sediment conductance as 1 m^2/d
    rbot = 0  # set location of river bottom (set to top of model, 0 m)

    # add rivers
    riv_spd={}
    for t in range(nper):
        rivers = []
        # for i in range(1,nrow-1):
        #     rivers.append([0, i, 0, stage_riv[t], cond, rbot])  # set left side as column of river cells

        # also add river cells for "horseshoe lake"
        for i in range(int(np.floor(nrow/3)-2), int(np.floor(nrow/3)+2)+1):
            for j in range(int(np.floor(ncol/2)-2), int(np.floor(ncol/2)+2)+1):
                rivers.append([0, i, j, lake_stage[t], cond, rbot])

        riv_spd[t] = rivers  # appropriate dictionary setup for flopy stress period data

    riv = flopy.modflow.ModflowRiv(m, stress_period_data=riv_spd, ipakcb=1)  # assign river package


    #%% ADD A WELL IN THE SW QUADRANT ========================================
    # add well
    wel_spd = {}
    # Create Single Well at center of domain with [lay, row, col, flux] list
    for t in range( nper ):
        pumping_rate = -100 + 10*t  # in m^3/d, negative for pumping/positive for injection
        well_1 = [0, 3*nrow/4, ncol/4, pumping_rate]
        wel_spd[t] = [well_1]
    # Create flopy wel object
    wel = flopy.modflow.ModflowWel(model=m, stress_period_data=wel_spd)


    #%% ESTABLISH INACTIVE CELLS ========================================
    #create ibound as array (1: active, 0: inactive, -1: constant head)
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    #assign NORTH and SOUTH of model domain as inactive
    ibound[:,0,:] = 0
    ibound[:,-1,:] = 0

    #designate starting heads as array of floats = 1.0 m
    strt = np.ones((nlay, nrow, ncol), dtype=np.float32)

    #create flopy "basic" object
    bas = flopy.modflow.ModflowBas(m, ibound=ibound, strt=strt)

    #%% DEFINE UPW PROPERTIES - UPW PACKAGE  ========================================
    #define horizontal and vertical hydraulic conductivity
    hk = np.ones((nlay,nrow,ncol), dtype=np.float32) * 100  # horizontal conductivity per cell = 1m/d
    vk = np.ones((nlay,nrow,ncol), dtype=np.float32) * 10  # vertical conductivity per cell = 1m/d

    #create flopy UPSTEAM WEIGHTING flow object
    upw = flopy.modflow.ModflowUpw(model=m,
                                   ipakcb=10, # writes the .cbc file!
                                   laytyp=1,
                                   hk=hk,
                                   layvka=1,
                                   vka=vk,
                                   sy=0.22,
                                   )

    #%% INITIATE RECHARGE RATE  ========================================
    rch_val = .005 #rch value per cell in m^3/d
    rch_array = np.zeros((nrow,ncol)) #rch array to hold values per cell

    #apply recharge to all cells except inactive and river cell boundary conditions
    rch_array[1:-1,1:-1] = rch_val
    rch_array[int(np.floor(nrow/3)-2):int(np.floor(nrow/3)+2),
              int(np.floor(ncol/2)-2):int(np.floor(ncol/2)+2) ]  = 0.0

    #create recharge object
    rch = flopy.modflow.ModflowRch(model=m,rech=rch_val) #EStL uses UZF package. need to update?

    # #%% CONSTANT HEAD DATA ========================================
    # chd_spd = {0: [[0, 21, 22, 6, 6]]}
    # chd = flopy.modflow.ModflowChd(m, stress_period_data=chd_spd)

    #%% WRITE THE MODEL FILES ========================================
    #write MODFLOW input files
    oc.reset_budgetunit() #<- reset the budget unit, because it is not consistently reset throughout the package calls.
    m.write_input()

    #%% CHECK THE MODEL BOUNDARY CONDITIONS - PRE-SIMULATION ========================================
    # inifitalize figure
    fig = plt.figure( figsize=(8,6) )
    #use flopy to check location of ibound, and rivers
    modelmap  = flopy.plot.PlotMapView(model=m, layer=0)
    grid      = modelmap.plot_grid()
    ib        = modelmap.plot_ibound()
    riv_cell  = modelmap.plot_bc(ftype='RIV',color='cyan')
    well_cell = modelmap.plot_bc(ftype='WEL',color='red')
    # chd_cell  = modelmap.plot_bc(ftype='CHD',color='yellow')
    #add labels and legend
    plt.xlabel('Lx (m)',fontsize = 14)
    plt.ylabel('Ly (m)',fontsize = 14)
    plt.title('Boundary Conditions', fontsize = 15, fontweight = 'bold')
    plt.legend(handles=[mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                        mp.patches.Patch(color='black',label='Inactive Cell',ec='black'),
                        mp.patches.Patch(color='cyan',label='Rivers',ec='black'),
                        mp.patches.Patch(color='red', label='Well',ec='black'),
                        mp.patches.Patch(color='yellow', label='Constant Head',ec='black'),
                        ], bbox_to_anchor=(1.01,0.5), loc='center left')
    plt.tight_layout()
    plt.show()


    #%% RUN THE MODEL ========================================
    import time
    start_time = time.time()
    # Run the model
    success, mfoutput = m.run_model(pause=False, report=True)
    run_time   = time.time() - start_time
    if not success:
        raise Exception('ISWS: MODFLOW did not terminate normally.')
    else:
        print(f'ISWS: MODFLOW model completed in {run_time} seconds.'+
              '\n'+'*****************************\n')

    # %% PLOT MODFLOW SIMULATION RESULTS ========================================
    from modules import headSurfacePlot, headContourPlot
    # identify the first time period
    totim = 1

    # obtain the head and flow budget objects
    headobj = flopy.utils.binaryfile.HeadFile(os.path.join(m._model_ws,model_name + '.hds'))
    budgobj = flopy.utils.binaryfile.CellBudgetFile(os.path.join(m._model_ws, model_name + '.cbc'))

    # get the data for the specific stress period
    head = headobj.get_data(totim=totim)
    figure = headSurfacePlot(m, Lx, Ly, head, title='Head Surface', fs=(8,8))

    # attempt flow figure
    frf = budgobj.get_data(text='flow right face', totim=totim)
    fff = budgobj.get_data(text='flow front face', totim=totim)
    fig2 = headContourPlot(m, head, frf, fff, title=' Flow arrows ', subtitle='Witness the greatness...')

    return m, run_time


