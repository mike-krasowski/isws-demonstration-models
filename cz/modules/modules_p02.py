"""
The purpose of this script is to house functions defined in the MODFLOW 6 examples repo in ex-prt-mp7-p02.py
"""


from pathlib import Path
from pprint import pformat

import flopy
import flopy.utils.binaryfile as bf
from flopy.plot.plotutil import bc_color_dict
# import git
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from flopy.mf6 import MFSimulation
from flopy.plot.styles import styles
from flopy.utils.gridgen import Gridgen
from flopy.utils.gridintersect import GridIntersect
from matplotlib.lines import Line2D
from modflow_devtools.misc import get_env, timed
from shapely.geometry import LineString, MultiPoint


def remove_COORDINATE_CHECK_METHOD(sim_ws):

    prt_path = sim_ws / 'prt'
    for item in os.listdir(prt_path):
        fn, fext = os.path.splitext(item)

        new_file = []
        if fext == '.prp':

            with open(prt_path / item, 'r') as fr:
                lines = fr.readlines()

                with open(prt_path / item, 'w') as fw:
                    for line in lines:
                        if "COORDINATE_CHECK_METHOD" not in line:
                            fw.write(line)

    print(
        """
                                    `;-.           ___,
               PRT ----->             `.`\\_....._/`.-"`
                                        \\         /       ,
                  ___________           /()    () \\     .' `-._
                 /    ow,    \\--/      ()   .    ()\\  /   _.'
                 \\  my self    /        \\  -'-     ,| '. <
                   -----------          ;.__      ,;|   > \\ 
                                       / ,    / ,   |.-'.-'
                                      (_/    (_/  ,;|.<`
                                        \\    ,      ;-`
                                         >   \\     /
                                        (_,-'`>  ./
                                              (_,/
        after jgs: https://www.asciiart.eu/art/ae81a937456fbb3d
        """
    )


def reverse_budgetfile(fpth, rev_fpth, tdis):
    f = bf.CellBudgetFile(fpth, tdis=tdis)
    f.reverse(rev_fpth)

def reverse_headfile(fpth, rev_fpth, tdis):
    f = bf.HeadFile(fpth, tdis=tdis)
    f.reverse(rev_fpth)

def get_mf6_pathlines(path):
    # load mf6 pathlines
    mf6pl = pd.read_csv(path)

    # index by particle group and particle ID
    mf6pl.set_index(["iprp", "irpt"], drop=False, inplace=True)

    # add column mapping particle group to subproblem name (1: A, 2: B)
    # mf6pl["subprob"] = mf6pl.apply(lambda row: "A" if row.iprp == 1 else "B", axis=1)

    # # add release time and termination time columns
    # mf6pl["t0"] = (
    #     mf6pl.groupby(level=["iprp", "irpt"])
    #     .apply(lambda x: x.t.min())
    #     .to_frame(name="t0")
    #     .t0
    # )
    # mf6pl["tt"] = (
    #     mf6pl.groupby(level=["iprp", "irpt"])
    #     .apply(lambda x: x.t.max())
    #     .to_frame(name="tt")
    #     .tt
    # )
    #
    # # add markercolor column, color-coding by layer for plots
    # mf6pl["mc"] = mf6pl.apply(
    #     lambda row: "green" if row.ilay == 1 else "yellow" if row.ilay == 2 else "red",
    #     axis=1,
    # )

    return mf6pl


def get_mp7_pathlines(path, gwf_model):
    # load mp7 pathlines, letting flopy determine capture areas
    mp7plf = flopy.utils.PathlineFile(path)
    mp7pl = pd.DataFrame(
        mp7plf.get_destination_pathline_data(
            list(range(gwf_model.modelgrid.nnodes)), to_recarray=True
        )
    )

    # index by particle group and particle ID
    mp7pl.set_index(["particlegroup", "sequencenumber"], drop=False, inplace=True)

    # convert indices to 1-based (flopy converts them to 0-based, but PRT uses 1-based, so do the same for consistency)
    kijnames = [
        "k",
        "node",
        "particleid",
        "particlegroup",
        "particleidloc",
        "sequencenumber",
    ]

    for n in kijnames:
        mp7pl[n] += 1

    # add column mapping particle group to subproblem name (1: A, 2: B)
    mp7pl["subprob"] = mp7pl.apply(
        lambda row: (
            "A" if row.particlegroup == 1 else "B" if row.particlegroup == 2 else pd.NA
        ),
        axis=1,
    )

    # add release time and termination time columns
    mp7pl["t0"] = (
        mp7pl.groupby(level=["particlegroup", "sequencenumber"])
        .apply(lambda x: x.time.min())
        .to_frame(name="t0")
        .t0
    )
    mp7pl["tt"] = (
        mp7pl.groupby(level=["particlegroup", "sequencenumber"])
        .apply(lambda x: x.time.max())
        .to_frame(name="tt")
        .tt
    )

    # add markercolor column, color-coding by layer for plots
    mp7pl["mc"] = mp7pl.apply(
        lambda row: "green" if row.k == 1 else "yellow" if row.k == 2 else "red", axis=1
    )

    return mp7pl


def get_mp7_timeseries(path, gwf_model):
    # load mp7 pathlines, letting flopy determine capture areas
    mp7tsf = flopy.utils.TimeseriesFile(path)
    mp7ts = pd.DataFrame(
        mp7tsf.get_destination_timeseries_data(list(range(gwf_model.modelgrid.nnodes)))
    )

    # index by particle group and particle ID
    mp7ts.set_index(["particlegroup", "particleid"], drop=False, inplace=True)

    # convert indices to 1-based (flopy converts them to 0-based, but PRT uses 1-based, so do the same for consistency)
    kijnames = ["k", "node", "particleid", "particlegroup", "particleidloc"]

    for n in kijnames:
        mp7ts[n] += 1

    # add column mapping particle group to subproblem name (1: A, 2: B)
    mp7ts["subprob"] = mp7ts.apply(
        lambda row: (
            "A" if row.particlegroup == 1 else "B" if row.particlegroup == 2 else pd.NA
        ),
        axis=1,
    )

    # add release time and termination time columns
    mp7ts["t0"] = (
        mp7ts.groupby(level=["particlegroup", "particleid"])
        .apply(lambda x: x.time.min())
        .to_frame(name="t0")
        .t0
    )
    mp7ts["tt"] = (
        mp7ts.groupby(level=["particlegroup", "particleid"])
        .apply(lambda x: x.time.max())
        .to_frame(name="tt")
        .tt
    )

    # add markercolor column, color-coding by layer for plots
    mp7ts["mc"] = mp7ts.apply(
        lambda row: "green" if row.k == 1 else "yellow" if row.k == 2 else "red", axis=1
    )

    return mp7ts


def get_mp7_endpoints(path, gwf_model):
    # load mp7 pathlines, letting flopy determine capture areas
    mp7epf = flopy.utils.EndpointFile(path)
    mp7ep = pd.DataFrame(
        mp7epf.get_destination_endpoint_data(list(range(gwf_model.modelgrid.nnodes)))
    )

    # index by particle group and particle ID
    mp7ep.set_index(["particlegroup", "particleid"], drop=False, inplace=True)

    # convert indices to 1-based (flopy converts them to 0-based, but PRT uses 1-based, so do the same for consistency)
    kijnames = ["k", "node", "particleid", "particlegroup", "particleidloc"]

    for n in kijnames:
        mp7ep[n] += 1

    # add column mapping particle group to subproblem name (1: A, 2: B)
    mp7ep["subprob"] = mp7ep.apply(
        lambda row: (
            "A" if row.particlegroup == 1 else "B" if row.particlegroup == 2 else pd.NA
        ),
        axis=1,
    )

    # add release time and termination time columns
    mp7ep["t0"] = (
        mp7ep.groupby(level=["particlegroup", "particleid"])
        .apply(lambda x: x.time.min())
        .to_frame(name="t0")
        .t0
    )
    mp7ep["tt"] = (
        mp7ep.groupby(level=["particlegroup", "particleid"])
        .apply(lambda x: x.time.max())
        .to_frame(name="tt")
        .tt
    )

    # add markercolor column, color-coding by layer for plots
    mp7ep["mc"] = mp7ep.apply(
        lambda row: "green" if row.k == 1 else "yellow" if row.k == 2 else "red", axis=1
    )

    return mp7ep


# -

# ### Plotting results
#
# Define functions to plot model results.


# +
# colormap for boundary locations
cmapbd = mpl.colors.ListedColormap(["r", "g"])

# time series point colors by layer
colors = ["green", "orange", "red"]

# figure sizes
figure_size_solo = (7, 7)
figure_size_compare = (7, 5)


# def plot_nodes_and_vertices(gwf, ax):
#     """
#     Plot cell nodes and vertices (and IDs) on a zoomed inset
#     """
#
#     ax.set_aspect("equal")
#
#     # set zoom area
#     xmin, xmax = 2050, 4800
#     ymin, ymax = 5200, 7550
#     ax.set_xlim([xmin, xmax])
#     ax.set_ylim([ymin, ymax])
#
#     # create map view plot
#     pmv = flopy.plot.PlotMapView(gwf, ax=ax)
#     pmv.plot_grid(edgecolor="black", alpha=0.25)
#     styles.heading(ax=ax, heading="Nodes and vertices (one-based)", fontsize=8)
#     ax.set_xlim([xmin, xmax])
#     ax.set_ylim([ymin, ymax])
#
#     # plot vertices
#     mg = gwf.modelgrid
#     verts = mg.verts
#     ax.plot(verts[:, 0], verts[:, 1], "bo", alpha=0.25, ms=2)
#     for i in range(ncpl):
#         x, y = verts[i, 0], verts[i, 1]
#         if xmin <= x <= xmax and ymin <= y <= ymax:
#             ax.annotate(str(i + 1), verts[i, :], color="b", alpha=0.5)
#
#     # plot nodes
#     xc, yc = mg.get_xcellcenters_for_layer(0), mg.get_ycellcenters_for_layer(0)
#     for i in range(ncpl):
#         x, y = xc[i], yc[i]
#         ax.plot(x, y, "o", color="grey", alpha=0.25, ms=2)
#         if xmin <= x <= xmax and ymin <= y <= ymax:
#             ax.annotate(str(i + 1), (x, y), color="grey", alpha=0.5)
#
#     # plot well
#     ax.plot(wel_coords[0][0], wel_coords[0][1], "ro")
#
#     # adjust left margin to compensate for inset
#     plt.subplots_adjust(left=0.45)
#
#     # create legend
#     ax.legend(
#         handles=[
#             Line2D(
#                 [0],
#                 [0],
#                 marker="o",
#                 color="w",
#                 label="Vertex",
#                 markerfacecolor="blue",
#                 markersize=10,
#             ),
#             Line2D(
#                 [0],
#                 [0],
#                 marker="o",
#                 color="w",
#                 label="Node",
#                 markerfacecolor="grey",
#                 markersize=10,
#             ),
#         ],
#         loc="upper left",
#     )


def plot_head(gwf, head, ibd, paths):
    with styles.USGSPlot():
        fig = plt.figure(figsize=figure_size_solo)
        fig.tight_layout()
        ax = fig.add_subplot(1, 1, 1, aspect="equal")
        # ax.set_xlim(0, Lx)
        # ax.set_ylim(0, Ly)
        ilay = 2
        cint = 0.25
        hmin = head[ilay, 0, :].min()
        hmax = head[ilay, 0, :].max()
        styles.heading(ax=ax, heading=f"Head, layer {ilay + 1}")
        mm = flopy.plot.PlotMapView(gwf, ax=ax, layer=ilay)
        mm.plot_bc("WEL", plotAll=True)
        mm.plot_bc("RIV", plotAll=True, color="teal")
        mm.plot_grid(alpha=0.25)

        # create inset
        # axins = ax.inset_axes([-0.8, 0.25, 0.7, 0.9])
        # # plot_nodes_and_vertices(gwf, axins)
        # ax.indicate_inset_zoom(axins)

        pc = mm.plot_array(head[:, 0, :], edgecolor="black", alpha=0.5)
        cb = plt.colorbar(pc, shrink=0.25, pad=0.1)
        cb.ax.set_xlabel(r"Head ($m$)")

        if ibd is not None:
            mm.plot_array(ibd, cmap=cmapbd, edgecolor="gray")

        # create legend
        ax.legend(
            handles=[
                mpl.patches.Patch(color="red", label="Well"),
                mpl.patches.Patch(color="teal", label="River"),
            ],
            loc="upper left",
        )

        levels = np.arange(np.floor(hmin), np.ceil(hmax) + cint, cint)
        cs = mm.contour_array(head[:, 0, :], colors="white", levels=levels)
        plt.clabel(cs, fmt="%.1f", colors="white", fontsize=11)

        # if plot_show:
        plt.show()
        # if plot_save:
        fig.savefig(paths['figs_path'] / f"{paths['sim_name']}-head")


def plot_points(ax, gwf, data):
    ax.set_aspect("equal")
    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(alpha=0.25)
    return ax.scatter(data["x"], data["y"],  s=3)  # color=data["mc"],


def plot_tracks(
    ax, gwf, title=None, ibd=None, pathlines=None, timeseries=None, endpoints=None
):
    ax.set_aspect("equal")
    # ax.set_xlim(0, Lx)
    # ax.set_ylim(0, Ly)
    if title is not None:
        ax.set_title(title, fontsize=12)

    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(alpha=0.25)

    if ibd is not None:
        mm.plot_array(ibd, cmap=cmapbd, edgecolor="gray")

    pts = []
    if pathlines is not None:
        pts.append(
            mm.plot_pathline(
                pathlines, layer="all", colors=["blue"], lw=0.5, ms=1, alpha=0.25
            )
        )
        if timeseries is None and endpoints is None:
            plot_points(ax, gwf, pathlines[pathlines.ireason != 1])
    if timeseries is not None:
        for k in range(0, gwf.modelgrid.nlay - 1):
            pts.append(mm.plot_timeseries(timeseries, layer=k, lw=0, ms=1))
            plot_points(ax, gwf, timeseries)
    if endpoints is not None:
        pts.append(
            mm.plot_endpoint(endpoints, direction="ending", colorbar=False, shrink=0.25)
        )
    return pts[0] if len(pts) == 1 else pts


def plot_pathlines_and_points(gwf, mf6pl, paths, title=None):
    fig, ax = plt.subplots(figsize=figure_size_compare)
    fig.tight_layout()

    plot_tracks(ax, gwf, paths, pathlines=mf6pl)

    axins = ax.inset_axes([-0.88, 0.2, 0.7, 0.9])
    plot_tracks(axins, gwf, paths, pathlines=mf6pl)
    ax.indicate_inset_zoom(axins)
    # xmin, xmax = 3000, 6300
    # ymin, ymax = 3500, 7000
    # axins.set_xlim([xmin, xmax])
    # axins.set_ylim([ymin, ymax])
    plt.subplots_adjust(left=0.55)

    ax.legend(
        title="EXPLANATION",
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                markersize=10,
                markerfacecolor="green",
                color="w",
                lw=4,
                label="Layer 1",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                markersize=10,
                markerfacecolor="gold",
                color="w",
                lw=4,
                label="Layer 2",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                markersize=10,
                markerfacecolor="red",
                color="w",
                lw=4,
                label="Layer 3",
            ),
        ],
    )

    if title is not None:
        styles.heading(ax, title)

    # if plot_show:
    plt.show()
    # if plot_save:
    fig.savefig(paths['figs_path'] / f"{paths['sim_name']}-paths")


def plot_2b(gwf, mf6endpoints, paths, title=None):
    fig, ax = plt.subplots(figsize=figure_size_compare)
    pts = plot_tracks(ax, gwf, endpoints=mf6endpoints)

    if title is not None:
        styles.heading(ax, title)

    cax = fig.add_axes([0.4, 0.12, 0.2, 0.01])
    cb = plt.colorbar(pts, cax=cax, orientation="horizontal", shrink=0.25)
    cb.set_label("Travel time (days)")

    plt.subplots_adjust(bottom=0.2)

    # if plot_show:
    plt.show()
    # if plot_save:
    fig.savefig(paths['figs_path'] / f"{paths['sim_name']}-endpts")


def plot_3d(gwf, pathlines, paths, sp=1, endpoints=None, title=None):
    import pyvista as pv
    from flopy.export.vtk import Vtk

    pv.set_plot_theme("document")
    axes = pv.Axes(show_actor=False, actor_scale=2.0, line_width=5)
    vert_exag = 10
    # pathlines = pathlines.to_records(index=False)
    vtk = Vtk(model=gwf, binary=False, vertical_exageration=vert_exag, smooth=False)
    vtk.add_model(gwf)
    vtk.add_pathline_points(pathlines)
    gwf_mesh, prt_mesh = vtk.to_pyvista()

    meshes = [gwf_mesh, prt_mesh]

    bc_meshes = {}
    for pckg_name in gwf.get_package_list():

        pckg = gwf.get_package(pckg_name)

        pckg_type = str(type(pckg))[-5:-2]

        if pckg.has_stress_period_data and (pckg_type != 'drn'):

            bc_meshes[pckg_type] = []

            print('ISWS: trblsht:', type(pckg))

            # when we assign a boundary condition once and let MODFLOW carry the condition forward without entering any
            # further information, this gets formatted weirdly. We need to count backward from the sp of interest to
            # find the most recent assignment that will be active.
            spda = pckg.stress_period_data.array
            spd = spda[sp]
            counter = 0
            while spda[sp - counter] is None:
                counter += 1
                spd = spda[sp-counter]
                if sp - counter < 0:
                    break

            if spd is not None:

                for bc_cell in spd:
                    bc_lrc = bc_cell[0]
                    print('ISWS: trblsht:', bc_lrc)

                    x_min = np.sum(gwf.modelgrid.delr[:bc_lrc[2]])
                    x_max = np.sum(gwf.modelgrid.delr[:bc_lrc[2] + 1])
                    y_min = np.sum(gwf.modelgrid.delc[bc_lrc[1] + 1:])
                    y_max = np.sum(gwf.modelgrid.delc[bc_lrc[1]:])
                    z_min = gwf.modelgrid.top_botm[bc_lrc[0]+1, bc_lrc[1], bc_lrc[2]] * vert_exag
                    z_max = gwf.modelgrid.top_botm[bc_lrc[0], bc_lrc[1], bc_lrc[2]] * vert_exag
                    bc_mesh = pv.Box(bounds=(x_min, x_max,
                                             y_min, y_max,
                                             z_min, z_max))

                    bc_meshes[pckg_type].append(bc_mesh)

                    meshes.append(bc_mesh)

    # now that we have all of the mesh objects in a list, loop through them and rotate to get our starting positions for
    # the figure
    for msh in meshes:
        msh.rotate_z(110, point=axes.origin, inplace=True)
        msh.rotate_y(-10, point=axes.origin, inplace=True)
        msh.rotate_x(10, point=axes.origin, inplace=True)

    plot_endpoints = False
    if endpoints is not None:
        plot_endpoints = True
        endpoints = endpoints.to_records(index=False)
        endpoints["z"] = endpoints["z"] * vert_exag
        eps_mesh = pv.PolyData(np.array(tuple(map(tuple, endpoints[["x", "y", "z"]]))))
        eps_mesh.rotate_z(-30, point=axes.origin, inplace=True)
        eps_mesh.rotate_y(-10, point=axes.origin, inplace=True)
        eps_mesh.rotate_x(10, point=axes.origin, inplace=True)

    def _plot(paths, screenshot=False):
        p = pv.Plotter(
            window_size=[500, 500],
            off_screen=screenshot,
            notebook=False if screenshot else None,
        )
        p.enable_anti_aliasing()
        if title is not None:
            p.add_title(title, font_size=7)
        p.add_mesh(gwf_mesh, opacity=0.025, style="wireframe")
        p.add_mesh(
            prt_mesh,
            # scalars="k" if "k" in prt_mesh.point_data else "ilay",
            # cmap=["green", "gold", "red"],
            point_size=4,
            line_width=3,
            render_points_as_spheres=True,
            render_lines_as_tubes=True,
            smooth_shading=True,
        )
        # p.remove_scalar_bar()
        if plot_endpoints:
            p.add_mesh(
                eps_mesh,
                # scalars=endpoints.k.ravel(),
                # cmap=["green", "gold", "red"],
                point_size=4,
            )
            # p.remove_scalar_bar()

        for key, bc_mesh in bc_meshes.items():

            # plot each of the cells using the standard color
            for bmsh in bc_mesh:
                p.add_mesh(bmsh, color=bc_color_dict[key.upper()], opacity=0.2)


        p.add_legend(
            labels=[("Layer 1", "green"), ("Layer 2", "gold"), ("Layer 3", "red")],
            bcolor="white",
            face="r",
            size=(0.15, 0.15),
        )

        p.camera.zoom(2)
        p.show()
        if screenshot:
            p.screenshot(paths['figs_path'] / f"{paths['sim_name']}-paths-3d.png")

    # if plot_show:
    _plot(paths)
    # if plot_save:
    _plot(paths, screenshot=True)


# def load_head():
#     head_file = flopy.utils.HeadFile(gwf_ws / (gwf_name + ".hds"))
#     return head_file.get_data()


# def get_ibound():
#     ibd = np.zeros((ncpl), dtype=int)
#     ibd[np.array(welcells)] = 1
#     ibd[np.array(rivcells)] = 2
#     return np.ma.masked_equal(ibd, 0)


def plot_all(gwf, paths):

    # ibound = get_ibound()
    mf6pl = get_mf6_pathlines(paths['prt_ws'] / paths['trackcsvfile_prt'])
    mf6ep = mf6pl[mf6pl.ireason == 3]  # termination event
    mp7pl = get_mp7_pathlines(paths['mp7_ws'] / f"{paths['mp7_name']}.mppth", gwf)
    mp7ts = get_mp7_timeseries(paths['mp7_ws'] / f"{paths['mp7_name']}.timeseries", gwf)
    # mp7ep = get_mp7_endpoints(paths['mp7_ws'] / f"{paths['mp7_name']}.mpend", gwf)

    # plot_head(gwf_model, load_head(), ibound)
    plot_pathlines_and_points(
        gwf,
        mf6pl=mf6pl,  # [mf6pl.subprob == "A"]
        paths=paths,
        title="Pathlines and points, 1000-day\ntime interval, colored by layer",
        # ibd=ibound,
    )
    plot_3d(
        gwf,
        pathlines=mf6pl,  # [mp7pl.subprob == "B"]
        paths=paths,
        endpoints=mf6ep,
        title="Pathlines and\nrecharge points,\ncolored by layer",
    )
    plot_2b(
        gwf,
        mf6endpoints=mf6ep,  # [mf6ep.subprob == "B"]
        paths=paths,
        title="Recharge points, colored by travel time",
        # ibd=ibound,
    )

