"""
The purpose of this script is to house a basic implementation of a PRT simulation in MODFLOW 6
"""

import flopy
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
import shutil

from modules.modules import *

# there are a lot of names, input files, and output files to keep track of in this workflow so it's best we define
# these early so they can be used throughout the workflow
sim_name = "0_basic"

gwf_name = sim_name + "-gwf"
prt_name = sim_name + "-prt"
mp7_name = sim_name + "-mp7"

# define an outputs folder
savepath = Path('outputs')

# if it exists, remove it and all of its contents
if os.path.exists(savepath):
    shutil.rmtree(savepath)

# make the outputs folder
os.makedirs(savepath)

sim_ws = savepath / sim_name
gwf_ws = sim_ws / "gwf"
prt_ws = sim_ws / "prt"
mp7_ws = sim_ws / "mp7"
figs_path = sim_ws / "figures"
gwf_ws.mkdir(exist_ok=True, parents=True)
prt_ws.mkdir(exist_ok=True, parents=True)
mp7_ws.mkdir(exist_ok=True, parents=True)
figs_path.mkdir(exist_ok=True, parents=True)

# Define output file names
gridfile = f"{gwf_name}.dis.grb"
headfile = f"{gwf_name}.hds"
budgetfile = f"{gwf_name}.cbb"
headfile_bkwd = f"{gwf_name}_bkwd.hds"
budgetfile_bkwd = f"{gwf_name}_bkwd.cbb"
budgetfile_prt = f"{prt_name}.cbb"
trackfile_prt = f"{prt_name}.trk"
trackhdrfile_prt = f"{prt_name}.trk.hdr"
trackcsvfile_prt = f"{prt_name}.trk.csv"
pathlinefile_mp7 = f"{mp7_name}.mppth"
endpointfile_mp7 = f"{mp7_name}.mpend"
timeseriesfile_mp7 = f"{mp7_name}.timeseries"

# MODFLOW 6 SIM - GWF ==================================================================================================
sim_gwf = flopy.mf6.MFSimulation(gwf_name, exe_name='../bin/win/mf6', sim_ws=gwf_ws)

# perlen, nstp, tsmult
perioddata = [(100000, 1, 1), ((5*365)+1, 1, 1)]

flopy.mf6.modflow.mftdis.ModflowTdis(
    sim_gwf,
    pname="tdis",
    nper=len(perioddata),
    perioddata=perioddata,
)

gwf = flopy.mf6.ModflowGwf(
    simulation=sim_gwf,
    modelname=gwf_name,
)

# iterative model solver (ims) package
ims = flopy.mf6.modflow.mfims.ModflowIms(
    sim_gwf,
    pname="ims",
    complexity="SIMPLE",
)
sim_gwf.register_solution_package(ims, [gwf.name])

nlay = 2
nrow = 10
ncol = 10
delr = 10
delc = 10
top = 100
botm = [50, 0]
elevs = [top] + botm

flopy.mf6.ModflowGwfdis(
    gwf,
    pname="dis",
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    delr=delr,
    delc=delc,
    top=top,
    botm=botm,
)

flopy.mf6.modflow.mfgwfic.ModflowGwfic(
    gwf,
    strt=top
)

kh = 1
kv = 0.1

flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
    gwf,
    # pname="npf",
    icelltype=1,
    k=kh,
    k33=kv,
    save_flows=True,
    save_specific_discharge=True,
    save_saturation=True,
)

sto = flopy.mf6.modflow.mfgwfsto.ModflowGwfsto(
    gwf,
    save_flows=True,
    iconvert=1,
    ss=0.0001,
    sy=0.1,
    steady_state={0: True},
    # transient={1: True},
)

flopy.mf6.ModflowGwfrcha(
    gwf,
    recharge=0.0001)

flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data={0:[[(0, 0, ccc), elevs[0]-20] for ccc in range(ncol)]+[[(0, nrow-1, ccc), elevs[1]+5] for ccc in range(ncol)]}
)

wel_spd = {1: [[(gwf.modelgrid.nlay - 1, int(gwf.modelgrid.nrow / 2), int(gwf.modelgrid.ncol / 2)), 100],
               [(gwf.modelgrid.nlay - 1, int(gwf.modelgrid.nrow / 3), int(gwf.modelgrid.ncol / 3)), 100]]}

flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(
    gwf,
    # maxbound=2,
    stress_period_data=wel_spd,
)

flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
    gwf,
    pname="oc",
    saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    head_filerecord=[f"{gwf.name}.hds"],
    budget_filerecord=[f"{gwf.name}.cbb"],
)

# MODFLOW 6 SIM - PRT ==================================================================================================
# simulation
sim_prt = flopy.mf6.MFSimulation(
    sim_name=prt_name, exe_name="../bin/win/mf6", version="mf6", sim_ws=prt_ws,
)

# temporal discretization
tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
    sim_prt,
    pname="tdis",
    # time_units=time_units,
    nper=len(perioddata),
    perioddata=perioddata,
)

# Instantiate the MODFLOW 6 prt model
prt = flopy.mf6.ModflowPrt(
    sim_prt, modelname=prt_name,
)

# Instantiate the MODFLOW 6 prt discretization package
flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
    prt,
    pname="dis",
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    length_units="FEET",
    delr=delr,
    delc=delc,
    top=top,
    botm=botm,
)

# Instantiate the MODFLOW 6 prt model input package.
flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=np.unique(gwf.sto.sy.array)[0] * 1.01)

well_lrc = wel_spd[1][0][0]
rel_minl = well_lrc[0]
rel_maxl = well_lrc[0]
rel_minr = well_lrc[1]
rel_maxr = well_lrc[1]
rel_minc = well_lrc[2]
rel_maxc = well_lrc[2]
celldata = flopy.modpath.CellDataType(
    drape=0,
    rowcelldivisions=5,
    columncelldivisions=5,
    layercelldivisions=10,
)
lrcregions = [[rel_minl, rel_minr, rel_minc, rel_maxl, rel_maxr, rel_maxc]]
lrcpd = flopy.modpath.LRCParticleData(
    subdivisiondata=[celldata],
    lrcregions=[lrcregions],
)
pg = flopy.modpath.ParticleGroupLRCTemplate(
    particlegroupname="PG1",
    particledata=lrcpd,
    # filename=f"{mp7_name}.pg1.sloc",
    # releasedata=(10, 0, 20),
)
pgs = [pg]
defaultiface = {"RECHARGE": 6, "ET": 6}
# Convert MODPATH 7 particle configuration to format expected by PRP.
release_points = list(lrcpd.to_prp(prt.modelgrid, localz=True))

# Instantiate the MODFLOW 6 prt particle release point (prp) package
pd = {0: ["FIRST"], 1: []}
flopy.mf6.ModflowPrtprp(
    prt,
    # pname='5yr',
    # filename=prpfilename,
    nreleasepts=len(release_points),
    packagedata=release_points,
    # perioddata=pd,
    exit_solve_tolerance=1e-5,
    extend_tracking=True,
)

# Instantiate the MODFLOW 6 prt output control package
track_times = list(range(90000, 150001, 2000))
flopy.mf6.ModflowPrtoc(
    prt,
    pname="oc",
    budget_filerecord=[budgetfile_prt],
    trackcsv_filerecord=[trackcsvfile_prt],
    track_release=True,
    track_terminate=True,
    track_usertime=True,
    ntracktimes=len(track_times),
    tracktimes=[(t,) for t in track_times],
    saverecord=[("BUDGET", "ALL")],
)

# Instantiate the MODFLOW 6 prt flow model interface
flopy.mf6.ModflowPrtfmi(
    prt,
    packagedata=[
        ("GWFGRID", Path(f"../{gwf_ws.name}/{gridfile}")),
        ("GWFHEAD", Path(f"../{gwf_ws.name}/{headfile}")),
        ("GWFBUDGET", Path(f"../{gwf_ws.name}/{budgetfile}")),
    ],
)

# Create an explicit model solution (EMS) for the MODFLOW 6 prt model
ems = flopy.mf6.ModflowEms(
    sim_prt,
    pname="ems",
    filename=f"{prt_name}.ems",
)
sim_prt.register_solution_package(ems, [sim_prt.name])

# MODPATH 7 ============================================================================================================
print("Building MODPATH 7 model...")

mp = flopy.modpath.Modpath7(
    modelname=mp7_name,
    flowmodel=gwf,
    exe_name="../bin/win/mpath7",
    model_ws=mp7_ws,
)
mpbas = flopy.modpath.Modpath7Bas(mp, porosity=np.unique(gwf.sto.sy.array)[0] * 1.01, defaultiface=defaultiface)
mpsim = flopy.modpath.Modpath7Sim(
    mp,
    pathlinefilename=pathlinefile_mp7,
    endpointfilename=endpointfile_mp7,
    timeseriesfilename=timeseriesfile_mp7,
    simulationtype="combined",
    trackingdirection="backward",
    weaksinkoption="pass_through",
    weaksourceoption="pass_through",
    budgetoutputoption="summary",
    # referencetime=[0, 0, 0.9],
    # timepointdata=[30, 2000.0],
    zonedataoption="on",
    # zones=izone,
    particlegroups=pgs,
)


# let's write all of our inputs
print('ISWS: writing the inputs for the *  gwf  * simulation')
sim_gwf.write_simulation(silent=False)
print('ISWS: writing the inputs for the *  prt  * simulation')
sim_prt.write_simulation(silent=False)
print('ISWS: writing the inputs for the *  mp7  * simulation')
mp.write_input()

# and now we're ready to run!
print('ISWS: running the *  gwf  * simulation')
success_gwf, _ = sim_gwf.run_simulation(silent=False, report=True)

# since we are running a "backwards" particle tracking scheme, reverse some important outputs
reverse_budgetfile(gwf_ws / budgetfile, gwf_ws / budgetfile_bkwd, sim_gwf.tdis)
reverse_headfile(gwf_ws / headfile, gwf_ws / headfile_bkwd, sim_gwf.tdis)

print('ISWS: running the *  prt  * simulation')
try:
    success_prt, _ = sim_prt.run_simulation(silent=False, report=True)
    assert success_prt
except:
    # PRT hurts itself in its confusion!
    remove_COORDINATE_CHECK_METHOD(sim_ws)
    success_prt, _ = sim_prt.run_simulation(silent=False, report=True)

print('ISWS: running the *  mp7  * simulation')
success, buff = mp.run_model(silent=False, report=True)


plot_all(
    gwf,
    paths=dict(
        sim_name = "0_basic",
        gwf_name = sim_name + "-gwf",
        prt_name = sim_name + "-prt",
        mp7_name = sim_name + "-mp7",
        savepath = Path('outputs'),
        sim_ws = savepath / sim_name,
        gwf_ws = sim_ws / "gwf",
        prt_ws = sim_ws / "prt",
        mp7_ws = sim_ws / "mp7",
        figs_path = sim_ws / "figures",
        gridfile = f"{gwf_name}.dis.grb",
        headfile = f"{gwf_name}.hds",
        budgetfile = f"{gwf_name}.cbb",
        headfile_bkwd = f"{gwf_name}_bkwd.hds",
        budgetfile_bkwd = f"{gwf_name}_bkwd.cbb",
        budgetfile_prt = f"{prt_name}.cbb",
        trackfile_prt = f"{prt_name}.trk",
        trackhdrfile_prt = f"{prt_name}.trk.hdr",
        trackcsvfile_prt = f"{prt_name}.trk.csv",
        pathlinefile_mp7 = f"{mp7_name}.mppth",
        endpointfile_mp7 = f"{mp7_name}.mpend",
        timeseriesfile_mp7 = f"{mp7_name}.timeseries",
    )
)
