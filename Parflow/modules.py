"""
The purpose of this script is to house functions for the layered groundwater model

author: ISWS, krasows2
July 16th, 2025
"""

# imports
import concurrent.futures
import numpy as np
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import os
import random
import shapely
from shapely import Polygon
import shutil
import time

topo_colors = ['lightgray', 'moccasin', 'tan', 'black', 'sandybrown', 'cornflowerblue']
tcmaps = {}
for tcolor in np.unique(topo_colors):
    name = 'alpha_to_' + tcolor
    tcmaps[name] = LinearSegmentedColormap.from_list(
        name,
        [(220/255, 0.0, 0.0, 0.0), colors.to_rgba(tcolor)],
        N=2)

colors = [(220/255, 0.0, 0.0, 0.0), (220/255, 0.0, 0.0, 1.0)]
atr_cmap = LinearSegmentedColormap.from_list('AlphaToRed', colors, N=100)

def find_cells_within_polygon( polygon, gridx, gridy ):
    """
    Find cells within the polygon. (adapted from ESL_MODFLOW_FloPy repo)

    Parameters
    ----------
    polygon : shapely.geometry
        Shapely.geometry polygon from a shapefile. A single geometry from a
        shapefile read in as a geopandas dataframe.
    gridx : numpy.ndarray
        Meshgrid of x-coordinates.
    gridy : numpy.ndarray
        Meshgrid of y-coordinates.

    Returns
    -------
    boolean array
        An array of boolean values designating whether a cell is within the
        provided polygon.

    """
    # create a boolean array of 'False' values (i.e., 0)
    pts_in = np.zeros( gridx.shape, dtype='bool')

    # get bounds of polygon <- [west, south, east, North]
    bbox = polygon.bounds
    # store the index of the array location where gridx, gridy are in bounds
    bbox_idx = np.where( (gridx >= bbox[0]) & (gridx <= bbox[2]) & \
                         (gridy >= bbox[1]) & (gridy <= bbox[3]) )

    for rowi, coli in zip( bbox_idx[0], bbox_idx[1] ):
        # check if polygon contains point -> boolean
        bool_val = polygon.contains(
                shapely.geometry.Point(gridx[rowi,coli],
                                       gridy[rowi,coli]) )
        if bool_val: # if true:
            # store True value
            pts_in[rowi, coli] = bool_val

    # ensuring points are found.
    if np.any( pts_in ):
        return pts_in
    else:
        return np.array([])
        print('Function did not find any cells within polygon.')


def determine_inside_indices(mf, polygons, axis=0, buffer=None):
    """
    The purpose of this function is to employ the find_cells_within_polygon() function to establish a field/column in
    a geopandas dataframe called "in_idx" which is a boolean array of which cells are within a polygon. The buffer
    parameter allows the user to expand the polygon in steps of the value supplied. This guarantees at least one model
    cell will be identified in association with each polygon. (adapted from ESL_MODFLOW_FloPy repo)

    :param mf: flopy.modflow.mf.Modflow
        a flopy model object which is used to get the grid coordinates
    :param polygons: geopandas.DataFrame
        geopandas dataframe containing polygons for which we want the associated model cells
    :param buffer: int, float, or None
        if int or float supplied, a buffer will be applied to a polygon if no cells are identified associated with it.
        if None is supplied (default) when no cells are identified, the script will simply move on.
    :return polygons: list or iterable
        list or iterable of shapely geometries for which inside modelgrid indices will be determined
    """

    # create a new field to store the original geometry
    polygons_og = polygons.copy()

    # create coordinate arrays of the desired shape

    list_of_in_idx = []
    for idx, polygon in enumerate(polygons):

        in_idx = np.array([])
        # across the layer plane (x-y)
        if axis == 0:

            # if buffer is supplied, expand the polygon to find cells if none are found in the original geometry
            if buffer is not None:
                while not in_idx.any():
                    in_idx = find_cells_within_polygon(polygon,
                                                       mf.modelgrid.xcellcenters,
                                                       mf.modelgrid.ycellcenters)
                    if not in_idx.any():
                        polygon = polygon.buffer(buffer)
            # else, use the standard function w/o buffering
            else:
                in_idx = find_cells_within_polygon(polygon,
                                                   mf.modelgrid.xcellcenters,
                                                   mf.modelgrid.ycellcenters)

        # across the row plane (x-z)
        elif axis == 1:

            # get the full shape of the x cell centers
            xcc = np.empty(mf.modelgrid.zcellcenters.shape)
            for idx in range(mf.modelgrid.zcellcenters.shape[0]):
                xcc[idx] = mf.modelgrid.xcellcenters.copy()

            # if buffer is supplied, expand the polygon to find cells if none are found in the original geometry
            if buffer is not None:
                while not in_idx.any():
                    in_idx = find_cells_within_polygon(polygon,
                                                       np.squeeze(xcc),
                                                       np.squeeze(mf.modelgrid.zcellcenters))
                    if not in_idx.any():
                        polygon = polygon.buffer(buffer)

            # else, use the standard function w/o buffering
            else:
                in_idx = find_cells_within_polygon(polygon,
                                                   np.squeeze(xcc),
                                                   np.squeeze(mf.modelgrid.zcellcenters))

        # across the column plane (y-z)
        elif axis == 2:

            # get the full shape of the y cell centers
            ycc = np.empty(mf.modelgrid.zcellcenters.shape)
            for idx in range(mf.modelgrid.zcellcenters.shape[0]):
                ycc[idx] = mf.modelgrid.ycellcenters.copy()

            # if buffer is supplied, expand the polygon to find cells if none are found in the original geometry
            if buffer is not None:
                while not in_idx.any():
                    in_idx = find_cells_within_polygon(polygon,
                                                       np.squeeze(ycc),
                                                       np.squeeze(mf.modelgrid.zcellcenters))
                    if not in_idx.any():
                        polygon = polygon.buffer(buffer)

            # else, use the standard function w/o buffering
            else:
                in_idx = find_cells_within_polygon(polygon,
                                                   np.squeeze(ycc),
                                                   np.squeeze(mf.modelgrid.zcellcenters))

        else:
            raise Exception('ISWS: Exception: axis argument not recognized, please enter 0, 1, or 2')

        # append result into list
        list_of_in_idx.append(in_idx)

    return list_of_in_idx

def get_nodes(gwf, locs):
    nodes = []
    for k, i, j in locs:
        nodes.append(k * gwf.modelgrid.nrow * gwf.modelgrid.ncol + i * gwf.modelgrid.ncol + j)
    return nodes

def run_the_models(sname, wiggle):

    # why do I have to do this here for pandas, but not for any other package.... they're all imported up top!
    import pandas as pd

    # import the geometries of the model (made by hand measurements of the table-top model
    topodata = pd.read_excel('./inputs/topology.xlsx')

    # create shapely polygons from measurements
    topopoly = {}
    for id in np.unique(topodata.id):
        feature = topodata[topodata.id == id]

        topopoly[id] = Polygon([(xxx, zzz) for xxx, zzz in zip(feature.x_coord, feature.z_coord)])

    # import the well information
    well_info = pd.read_csv('./inputs/well_info.csv')

    # build the flopy models represented by each scenario.
    #   The units of the model will be: length: cm, time: seconds, mass: milligrams

    print('ISWS: starting model for scenario:', sname)

    spds_path = './inputs/scenarios_short3.xlsx'
    exfi = pd.ExcelFile(spds_path)
    spd_schedule = exfi.parse(sname)
    exfi.close()

    perlen = np.array(spd_schedule.end - spd_schedule.start)


    # define a model workspace
    sim_ws = './outputs/{}'.format(sname)
    if not os.path.exists(sname):
        os.makedirs(sname)

    gwfname = 'gwf_' + sname
    
    # this is the main modflow model object. Many of the other objects interact
    sim = flopy.mf6.MFSimulation(
        sim_name='sim_' + sname,
        sim_ws=sim_ws,
        version="mf6",
        exe_name='../bin/win/mf6.exe'
    )

    flopy.mf6.ModflowTdis(
        sim,
        nper=spd_schedule.shape[0],
        perioddata=[[ddd, 1.0, 1.5] for ddd in perlen],
        time_units="SECONDS",
    )

    ### ----------------------------- specify solver -----------------------------

    # copied from: https://modflow6-examples.readthedocs.io/en/latest/_notebooks/ex-gwf-sagehen.html
    nouter, ninner = 10000, 1000
    hclose, rclose, relax = 1e-1, 1e-1, 0.97  # 3e-2, 3e-2, 0.97

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        # model_nam_file=f"{'gwf_' + sname}.nam",
        save_flows=True,
        newtonoptions="NEWTON UNDER_RELAXATION",
    )

    # preparing inputs for the geometry of the model (to be used by dis and npf)
    #   (inputs in the same order as topopoly dict)

    # dimensions of the model
    l_z = 29.8  # cm
    l_x = 51.2  # cm
    l_y = 2.54  # cm

    nlay = 60
    nrow = 1
    ncol = 100

    #  INITIALIZE DISCRETIZATION OBJECT
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=l_x / ncol,
        delc=l_y,
        top=l_z,
        botm=np.linspace(l_z - (l_z / nlay), 0, nlay),
        # idomain=1
    )

    in_idx_list = determine_inside_indices(gwf, list(topopoly.values()), axis=1, buffer=None)

    idomain = (in_idx_list[-3] * -1) + 1
    idomain[in_idx_list[-4]] = 0

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=l_x / ncol,
        delc=l_y,
        top=l_z,
        botm=np.linspace(l_z - (l_z / nlay), 0, nlay),
        idomain=idomain
    )

    #  SETTING INITIAL CONDITIONS
    flopy.mf6.ModflowGwfic(
        gwf,
        strt=l_z/1.5  # model top, doesn't have to be, though.
    )

    # INITIALIZE NODE PROPERTIES FILE
    hk_array = np.ones((gwf.modelgrid.nlay, gwf.modelgrid.nrow, gwf.modelgrid.ncol))
    for key, in_idx in zip(topopoly.keys(), in_idx_list):

        idx = np.where(in_idx)

        if key in ['aquitard', 'aquiclude_1', 'aquiclude_2', 'clay_layer_1', 'clay_layer_2']:

            hk_array[idx[0], :, idx[1]] = hk_array[idx[0], :, idx[1]] * 0.0001 / 600

        else:

            # the approximate speed of transport in a youtube video I saw (~20 cm in 6 mins)
            hk_array[idx[0], :, idx[1]] = hk_array[idx[0], :, idx[1]] * 0.5

    flopy.mf6.ModflowGwfnpf(
        gwf,
        # cvoptions="perched",
        # perched=True,
        icelltype=1,  # unconfined?
        k=hk_array,
        k33overk=True,
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

    # coordinate CHD and DRN locations at the top of the model
    riv_spd = {0: []}
    drn_spd = {0: []}
    chd_spd = {0: []}
    for ccc in range(gwf.modelgrid.ncol):

        for lll in range(idomain.shape[0]):
            if idomain[lll, ccc] == 1:
                break

        if ccc < 11:
            riv_spd[0].append([(lll, 0, ccc),
                               gwf.modelgrid.zcellcenters[20, 0, ccc] + 0.1,
                               2,
                               gwf.modelgrid.zcellcenters[lll, 0, ccc]])
        elif ccc <= 25:
            # do nothing
            pass

        elif ccc < 74:
            drn_spd[0].append([(lll, 0, ccc), gwf.modelgrid.zcellcenters[lll, 0, ccc], 1e5])

        else:
            chd_spd[0].append([(lll, 0, ccc), gwf.modelgrid.zcellcenters[lll, 0, ccc]+wiggle])

    # RIVERS RIVERS RIVERS
    riv = flopy.mf6.ModflowGwfriv(
        gwf,
        stress_period_data=riv_spd
    )

    # CONSTANT HEAD CONSTANT HEAD CONSTANT HEAD
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_spd
    )

    # DRAINS DRAINS DRAINS
    flopy.mf6.ModflowGwfdrn(
        gwf,
        stress_period_data=drn_spd
    )

    # get the location info (layer and column) of the elements of well_info
    well_layer = []
    well_row = []
    well_column = []
    for widx in well_info.index:

        # find the layer
        dell = gwf.modelgrid.top[0, 0] - gwf.modelgrid.botm[0, 0, 0]
        lll = np.round((l_z - well_info.screen_top[widx]) / dell).astype(int)

        # find the row
        rrr = 0

        # find the column
        ccc = np.round(well_info.x_coord[widx] / gwf.modelgrid.delr)[0].astype(int)

        well_layer.append(lll)
        well_row.append(rrr)
        well_column.append(ccc)

    well_info['lay'] = well_layer
    well_info['row'] = well_row
    well_info['col'] = well_column

    # WELLS WELLS WELLS
    q_factor = 0.1
    injection_conc = 100
    wel_spd = {}
    for sp in spd_schedule.index:
        wel_spd[sp] = []

        for widx in well_info.index:

            # find corresponding chem column
            for skey in spd_schedule.keys():
                if ('chem' in skey) and ('well' in well_info.id[widx]):
                    if skey[-1] == well_info.id[widx][-1]:

                        # set the injection well concentration
                        wel_conc = injection_conc * spd_schedule[skey][sp]

                        # get the pumping rate
                        qqq = spd_schedule[well_info.id[widx]][sp] * q_factor

                        # append into the dictionary
                        wel_spd[sp].append(
                            [(well_info.lay[widx], well_info.row[widx], well_info.col[widx]), qqq, wel_conc])

                    else:

                        # set the injection well concentration
                        wel_conc = 0

                        # get the pumping rate
                        qqq = spd_schedule[well_info.id[widx]][sp] * q_factor

                        # append into the dictionary
                        wel_spd[sp].append(
                            [(well_info.lay[widx], well_info.row[widx], well_info.col[widx]), qqq, wel_conc])

    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=wel_spd,
        auxiliary="CONCENTRATION",
        save_flows=False,
        maxbound=well_info.shape[0] * 2,
        pname="WEL-1",
    )

    head_filerecord = f"{gwf.name}.hds"
    budget_filerecord = f"{gwf.name}.cbc"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    nouter = 10000
    ninner = 10000
    outer_dvclose = 1e-3,
    inner_dvclose = 1e-3,
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        complexity='complex',
        inner_maximum=ninner,
        outer_maximum=nouter,
        outer_dvclose=outer_dvclose,
        inner_dvclose=inner_dvclose,
    )

    sim.register_ims_package(imsgwf, [gwf.name])

    # gwf.header = (' ISWS: “A purpose of human life, no matter who is controlling it, ' +
    #               'is to love whoever is around to be loved.” - Vonnegut')

    # TRANSPORT TRANSPORT TRANSPORT

    gwtname = gwf.name.replace('gwf', 'gwt')

    # start the gwt model object
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )

    gwt.name_file.save_flows = True

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=gwf.modelgrid.nlay,
        nrow=gwf.modelgrid.nrow,
        ncol=gwf.modelgrid.ncol,
        delr=gwf.modelgrid.delr,
        delc=gwf.modelgrid.delc,
        top=gwf.modelgrid.top,
        botm=gwf.modelgrid.botm,
        idomain=gwf.dis.idomain.array,
        filename=f"{gwt.name}.dis",
    )

    # starting concentrations
    sconc = np.zeros((gwf.modelgrid.nlay, gwf.modelgrid.nrow, gwf.modelgrid.ncol), dtype=float)

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGwtic(
        gwt,
        strt=sconc,
        filename=f"{gwt.name}.ic"
    )

    # Instantiating MODFLOW 6 transport advection package
    # HMOC parameters in case they are invoked in the ADV package
    dceps = 1.0e-5  # HMOC
    nplane = 1  # HMOC
    npl = 0  # HMOC
    nph = 4  # HMOC
    npmin = 0  # HMOC
    npmax = 8  # HMOC
    nlsink = nplane  # HMOC
    npsink = nph  # HMOC

    # if adv_scheme_mixelm == 0:
    #     scheme = "UPSTREAM"
    # elif adv_scheme_mixelm == -1:
    #     scheme = "TVD"
    # else:
    #     raise Exception()

    scheme = "TVD"
    flopy.mf6.ModflowGwtadv(
        gwt,
        scheme=scheme,
        filename=f"{gwt.name}.adv"
    )

    # Instantiating MODFLOW 6 transport dispersion package
    dsp_dispersivity = 0.1
    dsp_dmcoef = 1e-6
    flopy.mf6.ModflowGwtdsp(
        gwt,
        diffc=dsp_dmcoef,
        alh=dsp_dispersivity,
        ath1=dsp_dispersivity * 0.1,
        atv=dsp_dispersivity * 0.01,
        pname="DSP-1",
        filename=f"{gwt.name}.dsp",
    )

    # Instantiating MODFLOW 6 transport mass storage package (formerly "reaction" package in MT3DMS)
    prsity = gwf.sto.sy.array
    # if mst_kd_retardation != 1.0:
    #     sorption = "linear"
    #     bulk_density = mst_kd_rhob
    #     kd = (mst_kd_retardation - 1.0) * prsity / mst_kd_rhob
    # else:  # global variable section
    #     sorption = None
    #     bulk_density = None
    #     kd = None
    # if mst_decay != 0.0:
    #     first_order_decay = True
    #     decay_arg = mst_decay
    # else:
    #     first_order_decay = False
    #     decay_arg = None

    sorption = None
    bulk_density = None
    kd = None
    first_order_decay = False
    decay_arg = None

    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=prsity,
        sorption=sorption,
        bulk_density=bulk_density,
        distcoef=kd,
        first_order_decay=first_order_decay,
        decay=decay_arg,
        decay_sorbed=decay_arg,
        filename=f"{gwt.name}.mst",
    )

    # this next bit only works if we are only introducing chem one cell at a time
    cnc_conc = injection_conc * 0.005
    cnc_spd = {sp:[] for sp in spd_schedule.index}

    linger_on = True
    for key in spd_schedule.keys():
        linger = False
        for sp in spd_schedule.index:
            if 'chem' in key:
                if (spd_schedule.loc[sp, key] == 1) or linger:
                    if linger_on:
                        linger = True
                    for widx in well_info.index:
                        if ('well' in well_info.id[widx]) and (key[-1] == well_info.id[widx][-1]):
                            cnc_spd[sp].append(
                                [(well_info.lay[widx], well_info.row[widx], well_info.col[widx]), cnc_conc]
                            )

    maxbound = 0
    for key, item in cnc_spd.items():
        if len(item) > maxbound:
            maxbound = len(item)

    cnc = flopy.mf6.ModflowGwtcnc(
        gwt,
        # boundnames=True,
        maxbound=maxbound+1,
        stress_period_data=cnc_spd,
        save_flows=False,
        pname="CNC-1",
        filename=f"{gwt.name}.cnc",
    )

    # initialize source sink mixing package
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=sourcerecarray,
        filename=f"{gwt.name}.ssm"
    )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwf.name,
        exgmnameb=gwt.name,
        filename=f"{gwf.name + '-' + gwt.name}.gwfgwt",
    )
    # pkgdata = [
    #     ("GWFHEAD", gwf.name + '.hds', None),
    #     ("GWFBUDGET", gwf.name + '.bud', None),
    # ]
    # flopy.mf6.ModflowGwtfmi(
    #     gwt,
    #     flow_imbalance_correction=True,
    #     packagedata=pkgdata
    # )

    # this is only needed if you want to use separate settings for the flow and transport models
    hclose, rclose, relax = 1e-2, 1e-1, 0.97  # 3e-2, 3e-2, 0.97

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity='complex',
        outer_dvclose=outer_dvclose,
        inner_dvclose=inner_dvclose,
        # under_relaxation="NONE",
        inner_maximum=ninner,
        outer_maximum=nouter,
        # inner_dvclose=hclose,
        # rcloserecord=rclose,
        # linear_acceleration="BICGSTAB",
        # scaling_method="NONE",
        # reordering_method="NONE",
        # relaxation_factor=relax,
        filename="{}.ims".format("gwtsolver"),
    )

    sim.register_ims_package(imsgwt, [gwt.name])

    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwt.name}.cbc",
        concentration_filerecord=f"{gwt.name}.ucn",
        concentrationprintrecord=[("COLUMNS", gwf.modelgrid.ncol, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    # mfsetup/flopy/MODFLOW 6 expects the model objects to have some attributes which the transport models do not have.
    #   --> We can fake them!
    # gwt.header = ' ISWS: “I admire anybody who finishes a work of art, no matter how awful it may be.” - Vonnegut'
    #
    # gwt.cfg = dict(mfsetup_options=dict(keep_original_arrays=True))
    #
    # gwt.name_file.save_flows = True

    # WRITE THE INPUTS
    sim.write_simulation()

    # RUN THE FLOW AND TRANSPORT MODELS
    success_mf, buff = sim.run_simulation(silent=False)

    # MODPATH MODPATH MODPATH
    mp_path = os.path.join(sim.sim_path, 'mp')
    if not os.path.exists(mp_path):
        os.makedirs(mp_path)

    mp_locs = []

    for lll in range(gwf.modelgrid.nlay):
        for rrr in range(gwf.modelgrid.nrow):
            for ccc in range(gwf.modelgrid.ncol):

                if idomain[lll, ccc] == 1:

                    if (lll % 5 == 0) and (ccc % 5 == 0) and np.random.choice([True, True, True, False, False]):

                        lllw = lll + np.random.choice([-2, -1, 0, 1, 2])  # [-2, -1, 0, 1, 2]
                        if lllw >= gwf.modelgrid.nlay:
                            lllw = gwf.modelgrid.nlay - 1
                        elif lllw < 0:
                            lllw = 0

                        rrrw = rrr + np.random.choice([-2, -1, 0, 1, 2])
                        if rrrw >= gwf.modelgrid.nrow:
                            rrrw = gwf.modelgrid.nrow - 1
                        elif rrrw < 0:
                            rrrw = 0

                        cccw = ccc + np.random.choice([-2, -1, 0, 1, 2])
                        if cccw >= gwf.modelgrid.ncol:
                            cccw = gwf.modelgrid.ncol - 1
                        elif cccw < 0:
                            cccw = 0

                        mp_locs.append((int(lllw), int(rrrw), int(cccw)))

    # for sp, entry in chd_spd.items():
    #     for ch in entry:
    #         mp_locs.append(ch[0])

    # mp_locs.append((20, 0, 95))
    # mp_locs.append((25, 0, 95))
    # mp_locs.append((30, 0, 95))
    # mp_locs.append((35, 0, 95))
    # mp_locs.append((40, 0, 95))
    # mp_locs.append((45, 0, 95))

    nodes = gwf.modelgrid.get_node(mp_locs)
    
    # create basic forward tracking modpath simulation
    mp = flopy.modpath.Modpath7(
        modelname=gwf.name.replace('gwf', 'mp'),
        flowmodel=gwf,
        model_ws=mp_path,
        exe_name='../bin/win/mpath7.exe',
    )

    flopy.modpath.Modpath7Bas(mp, porosity=gwf.sto.sy.array)  # defaultiface=defaultiface)

    cd = flopy.modpath.CellDataType(
        drape=0,  # particles added at top of cell (no drape), 1 means they fall to active cell
        rowcelldivisions=1,
        columncelldivisions=1,
        layercelldivisions=1,
    )

    pd = flopy.modpath.NodeParticleData(
        subdivisiondata=[cd],
        nodes=nodes
    )

    # release releasedata[0] times starting at releasedata[1] and do so every releasedata[2] stress periods
    pg = flopy.modpath.ParticleGroupNodeTemplate(
        particlegroupname='PG1',
        # filename=None,
        releasedata=[gwf.nper//10, 1, 10],
        particledata=pd)

    pgs = [pg]
    denominator = 5
    number_of_times_particles_are_introduced = np.floor(gwf.nper / denominator).astype(int)
    time_between_introduction = denominator * 1
    # release how many times, with how long between releases
    timepointdata = [number_of_times_particles_are_introduced,
                     [time_between_introduction] * number_of_times_particles_are_introduced]
    timepointdata = list(spd_schedule.loc[1:, 'start']) + [spd_schedule.loc[len(spd_schedule) - 1, 'end']]

    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="combined",
        trackingdirection="forward",
        weaksinkoption="pass_through",
        weaksourceoption="pass_through",
        budgetoutputoption="summary",
        # referencetime=[0, 0, 0.9],
        # timepointdata=[len(timepointdata), timepointdata],
        # zonedataoption="on",
        # zones=zone_maps,
        particlegroups=pgs,
    )

    # write modpath datasets
    mp.write_input()

    # run modpath
    success_mp, buff = mp.run_model(silent=True, report=True)
    for line in buff:
        print(line)

    success = False
    if success_mf and success_mp:
        success = True
    else:
        print('ISWS: success_mf is {} and success_mp is {}.'.format(success_mf, success_mp))

    return success, sim, nodes


def make_the_animation(sim, nodes):

    # animation building will go here and it will be built using the scd variable from above
    savepath = './outputs/animations'
    if not os.path.exists(savepath):
        os.makedirs(savepath)

    FFMpegWriter = animation.writers['ffmpeg']
    vid = FFMpegWriter(fps=4)

    names = list(sim._models.keys())

    gwf = sim.get_model(names[0])
    gwt = sim.get_model(names[1])

    sname = gwf.name[-2:]

    # import the geometries of the model (made by hand measurements of the table-top model
    topodata = pd.read_excel('./inputs/topology.xlsx')

    # create shapely polygons from measurements
    topopoly = {}
    for id in np.unique(topodata.id):
        feature = topodata[topodata.id == id]

        topopoly[id] = Polygon([(xxx, zzz) for xxx, zzz in zip(feature.x_coord, feature.z_coord)])

    in_idx_list = determine_inside_indices(gwf, list(topopoly.values()), axis=1, buffer=None)

    print('ISWS: gwf.name:', gwf.name)

    print('ISWS: starting animation for scenario:', sname)

    fig, ax = plt.subplots()

    spds_path = './inputs/scenarios_short3.xlsx'
    exfi = pd.ExcelFile(spds_path)
    spd_schedule = exfi.parse(sname)
    exfi.close()

    # load heads
    fname = os.path.join(gwf.model_ws[:-1], gwf.name + '.hds')
    head_obj = flopy.utils.binaryfile.HeadFile(fname)
    heads = head_obj.get_alldata()
    head_obj.close()

    # swap the inactive value with np.nan
    heads = np.where(heads>=1e30, np.nan, heads)

    # load concentrations
    concs = gwt.output.concentration().get_alldata()

    input_conc = 100
    concs = np.where(concs > input_conc, np.nan, concs)

    # load paths
    fpth = os.path.join(sim.sim_path, 'mp', f"mp_{sname}.mppth")
    pathline_file = flopy.utils.PathlineFile(fpth)
    pathlines = pathline_file.get_destination_pathline_data(dest_cells=nodes)  # to_recarray=True

    # raise Exception

    void_dtype = pathlines[0][0].dtype

    # pre-process the path lines
    pathlines_by_sp = {}
    pathlines_by_spv = {}
    for sp in spd_schedule.index:
        pathlines_by_sp[sp] = []
        pathlines_by_spv[sp] = []
        for pthl in pathlines:

            pathlines_by_sp[sp].append(([], [], []))
            pathlines_by_spv[sp].append(np.array([], dtype=void_dtype))

            # each point in the path is a np.void object
            for pvoid in pthl:

                if (pvoid['stressperiod'] == sp) or (pvoid['stressperiod'] == sp - 1):
                    pathlines_by_sp[sp][-1][0].append(pvoid['x'])
                    pathlines_by_sp[sp][-1][1].append(pvoid['y'])
                    pathlines_by_sp[sp][-1][2].append(pvoid['z'])

                    pathlines_by_spv[sp].append(pvoid)
                    pathlines_by_spv[sp][-1] = np.append(pathlines_by_spv[sp][-1], pvoid)

    # also accepts .gif format which is helpful for filling out the readme with examples
    with vid.saving(fig, os.path.join(savepath, '{}.mp4'.format(sname)), 600):

        for sp in range(gwf.nper):

            if np.remainder(sp, 5) == 0:
                print('ISWS: {}: Drawing stress period: {}'.format(sname, sp))

            # plot data
            xsec_r = flopy.plot.PlotCrossSection(model=gwf,
                                                 line={'Row': 0},
                                                 ax=ax)

            xsec_r.plot_grid(linewidths=0.1, color='black', zorder=10000)

            for pidx, (id, poly) in enumerate(topopoly.items()):

                lith_array = in_idx_list[pidx]

                if id in ['aquitard', 'aquiclude_1', 'aquiclude_2',
                          'clay_layer_1', 'clay_layer_2', 'fractured_bedrock']:
                    xsec_r.plot_array(lith_array, cmap=tcmaps['alpha_to_lightgray'])
                elif id in ['soil_zone', 'confined_artesian_aquifer']:
                    xsec_r.plot_array(lith_array, cmap=tcmaps['alpha_to_moccasin'])
                elif id in ['fine_sand']:
                    xsec_r.plot_array(lith_array, cmap=tcmaps['alpha_to_tan'])
                elif id in ['unconfined_aquifer']:
                    xsec_r.plot_array(lith_array, cmap=tcmaps['alpha_to_sandybrown'])
                elif id in ['inactive']:
                    xsec_r.plot_array(lith_array, cmap=tcmaps['alpha_to_black'])
                elif id in ['river_channel']:
                    xsec_r.plot_array(lith_array, cmap=tcmaps['alpha_to_cornflowerblue'])
                else:
                    print('ruh roh!')

            # plot wells
            well_info = pd.read_csv('./inputs/well_info.csv')

            for widx in well_info.index:
                wid = well_info.id[widx]
                if wid in list(spd_schedule.keys()):
                    if spd_schedule[wid][sp] != 0:
                        wcolor = (1, 0, 0, 0.75)
                    else:
                        wcolor = (1, 1, 1, 0.5)
                else:
                    wcolor = (1, 1, 1, 0.5)

                w_x1 = well_info.x_coord[widx]
                w_z1 = well_info.z_coord[widx]
                w_z2 = well_info.screen_top[widx]

                if not np.isnan(well_info.screen_bottom[widx]):
                    w_z2 = well_info.screen_bottom[widx]

                ax.plot([w_x1, w_x1],[w_z1, w_z2], color=wcolor, linewidth=2)

            xsec_r.plot_array(concs[sp], cmap=atr_cmap, zorder=5000)

            # plot pathlines
            if sp >= 1:
                for pthl in pathlines_by_sp[sp]:
                    if np.random.choice([True, True, True, True, False]):
                        ax.plot(pthl[0], pthl[2], color='cornflowerblue', lw=0.75, zorder=1000)
                        # counter += 1

            # plot potentiometric surface
            pot_surf = np.zeros(heads.shape[3])
            for ccc in range(heads.shape[3]):

                for lll in range(gwf.dis.idomain.array.shape[0]):
                    if gwf.dis.idomain.array[lll, 0, ccc] == 1:
                        break


                pot_surf[ccc] = np.nanmax(heads[sp, lll, 0, ccc])

            ax.plot(gwf.modelgrid.xcellcenters.flatten(), pot_surf, lw=1.25, color='cornflowerblue')

            # figure finagling
            ax.set_title('Scenario: {}, time: {}'.format(sname, spd_schedule.end[sp]))

            ax.axis('equal')
            ax.axis('off')

            fig.tight_layout()


            # fig.legend(handles=hdls, loc='upper left', bbox_to_anchor=(0.81, 0.62), title='Concentration, (ppb)')
            # set super title for figure
            # fig.suptitle('End of MF Stress Period {}'.format(sp))

            # grab the figure before closing it
            vid.grab_frame()
            # clear the axes?

            ax.clear()

        # finish making movie?
        vid.finish()
        plt.close()
        print('Animation complete.')
        print('Saved as:\n     >> {}'.format(savepath))

    return sim