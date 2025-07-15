"""
The purpose of this script is to build MODFLOW models to simulate the various scenarios of the table-top layered
groundwater model and produce animations of the results

author: ISWS, krasows2, zavelle
July 14th, 2025
"""

# imports
import flopy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# define an outputs folder
savepath = 'outputs'
if not os.path.exists(savepath):
    os.makedirs(savepath)

exfi = pd.ExcelFile('./inputs/scenarios.xlsx')
scenarios = exfi.sheet_names

# build the flopy models represented by each scenario
scd = {}
for scenario in scenarios:

    spd_schedule = exfi.parse(scenario)

    model_name = scenario

    # define a model workspace
    sim_ws = './outputs/{}'.format(scenario)
    if not os.path.exists(savepath):
        os.makedirs(savepath)

    # this is the main modflow model object. Many of the other objects interact
    sim = flopy.mf6.MFSimulation(
        sim_name=model_name,
        sim_ws=sim_ws,
        version="mf6",
        exe_name='../bin/win/mf6.exe'
    )

    flopy.mf6.ModflowTdis(
        sim,
        nper=spd_schedule.shape[0],
        perioddata=[[ddd, 1.0, 1.5] for ddd in spd_schedule.end - spd_schedule.start],
        time_units="SECONDS",
    )

    ### ----------------------------- specify solver -----------------------------

    # copied from: https://modflow6-examples.readthedocs.io/en/latest/_notebooks/ex-gwf-sagehen.html
    nouter, ninner = 300, 500
    hclose, rclose, relax = 3e-2, 3e-2, 0.97

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        # outer_dvclose=hclose,
        # outer_maximum=nouter,
        # under_relaxation="dbd",
        linear_acceleration="BICGSTAB",
        # under_relaxation_theta=0.7,
        # under_relaxation_kappa=0.08,
        # under_relaxation_gamma=0.05,
        # under_relaxation_momentum=0.0,
        # inner_dvclose=hclose,
        # rcloserecord="1000.0 strict",
        # inner_maximum=ninner,
        # relaxation_factor=relax,
        # number_orthogonalizations=2,
        # preconditioner_levels=8,
        # preconditioner_drop_tolerance=0.001,
        filename=f"{model_name}.ims",
    )
    # sim.register_ims_package(imsgwf, [model_name])

    model_nam_file = f"{model_name}.nam"

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=model_name,
        model_nam_file=model_nam_file,
        save_flows=True,
        newtonoptions="NEWTON UNDER_RELAXATION",
    )

    #  INITIALIZE DISCRETIZATION OBJECT
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=60,
        nrow=1,
        ncol=100,
        delr=1,
        delc=1,
        top=100,
        botm=np.linspace(100-(100/60), 0, 60)
    )

    #  SETTING INITIAL CONDITIONS
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=100
    )

    # INITIALIZE NODE PROPERTIES FILE
    hk_array = np.ones((gwf.modelgrid.nlay, gwf.modelgrid.ncol))
    hk_array = hk_array*(50/600)  # the approximate speed of transport in a youtube video I saw

    flopy.mf6.ModflowGwfnpf(
        gwf,
        # cvoptions="perched",
        # perched=True,
        icelltype=1,  # unconfined?
        k=hk_array,
        k33=1,
        save_specific_discharge=True,
    )

    # INITIALIZE STORAGE

    steady_state = {0: True}
    transient = {1: True}

    flopy.mf6.ModflowGwfsto(
        gwf,
        pname="sto",
        save_flows=True,
        iconvert=1,
        # ss=ss_array,
        # sy=sy_array,
        steady_state=steady_state,
        transient=transient,
    )

    # RECHARGE RECHARGE RECHARGE
    flopy.mf6.ModflowGwfrcha(
        gwf,
        recharge=1/100
    )

    # DRAINS DRAINS DRAINS
    drn_spd = {0:[]}
    for ccc in range(gwf.modelgrid.ncol):

        if ccc < 5:

            drn_spd[0].append([0, 0, ccc, 99, 1e5])

        else:

            drn_spd[0].append([0, 0, ccc, 100, 1e5])

    # initialize the drains object
    flopy.mf6.ModflowGwfdrn(
        gwf,
        stress_period_data=drn_spd
    )


    head_filerecord = f"{model_name}.hds"
    budget_filerecord = f"{model_name}.cbc"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # WRITE THE INPUTS
    sim.write_simulation()

    # RUN THE MODEL
    success, buff = sim.run_simulation(silent=False)
    if not success:
        raise Exception("MODFLOW 6 did not terminate normally.")

    scd[scenario] = sim


# animation building will go here and be built using the scd variable from above


