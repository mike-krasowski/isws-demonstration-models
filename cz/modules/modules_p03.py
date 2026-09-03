"""
The purpose of this script is to house modules in support of the main series of scripts in this folder starting with
0_basic.py
"""

import flopy
from flopy.export.vtk import Vtk
import flopy.utils.binaryfile as bf
from flopy.plot.styles import styles
from flopy.plot.plotutil import bc_color_dict
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import pandas as pd
import pyvista as pv

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

# +
# Pathline and starting point colors by destination
colordest = {"well": "red", "drain": "green", "river": "blue"}


def plot_pathlines(ax, gwf, data, **kwargs):
    ax.set_aspect("equal")
    pl = flopy.plot.PlotMapView(model=gwf, ax=ax)
    pl.plot_grid(lw=0.5, alpha=0.5)
    pl.plot_bc("WEL", plotAll=True)
    pl.plot_bc("RIV", plotAll=True)
    pl.plot_bc("DRN", plotAll=True, color="green")

    for dest in ["well", "drain", "river"]:
        label = None
        if "colordest" in kwargs:
            color = kwargs["colordest"][dest]
            label = "Captured by " + dest
        elif "color" in kwargs:
            color = kwargs["color"]
        else:
            color = "grey"
        d = data[data.dest == dest]
        pl.plot_pathline(
            d, layer="all", colors=[color], label=label, linewidth=0.5, alpha=0.5
        )


def plot_points(fig, ax, gwf, data, colorbar=True, **kwargs):
    ax.set_aspect("equal")
    mm = flopy.plot.PlotMapView(model=gwf, ax=ax)
    mm.plot_grid(lw=0.5, alpha=0.5)
    if "colordest" in kwargs:
        pts = []
        for dest in ["well", "river", "drain"]:
            color = kwargs["colordest"][dest]
            label = "Captured by " + dest
            pdata = data[data.dest == dest]
            pts.append(
                ax.scatter(pdata["x"], pdata["y"], s=3, color=color, label=label)
            )
        return pts
    else:
        pts = ax.scatter(
            data["x"],
            data["y"],
            s=3,
            c=data["t"] if "t" in data.dtypes.keys() else data["time"],
        )
        if colorbar:
            cax = fig.add_axes([0.2, 0.19, 0.6, 0.01])
            cb = plt.colorbar(pts, cax=cax, orientation="horizontal", shrink=0.25)
            cb.set_label("Travel time (days)")


def plot_head(gwf, paths, head):
    with styles.USGSPlot():
        fig, ax = plt.subplots(figsize=(7, 7))
        fig.tight_layout()
        ax.set_aspect("equal")
        ilay = -1
        cint = 0.25
        hmin = head[ilay, 0, :].min()
        hmax = head[ilay, 0, :].max()
        styles.heading(ax=ax, heading=f"Head, layer {ilay + 1!s}, time=0")
        mm = flopy.plot.PlotMapView(gwf, ax=ax, layer=ilay)
        mm.plot_grid(lw=0.5)
        mm.plot_bc("WEL", plotAll=True)
        # mm.plot_bc("CHD", plotAll=True)
        mm.plot_bc("RIV", plotAll=True)
        # mm.plot_bc("DRN", plotAll=True, color="green")

        pc = mm.plot_array(head[ilay, :, :], edgecolor="black", alpha=0.25)
        cb = plt.colorbar(pc, shrink=0.25, pad=0.1)
        cb.ax.set_xlabel(r"Head ($ft$)")

        levels = np.arange(np.floor(hmin), np.ceil(hmax) + cint, cint)
        cs = mm.contour_array(head[ilay, :, :], colors="white", levels=levels)
        plt.clabel(cs, fmt="%.1f", colors="white", fontsize=11)

        ax.legend(
            handles=[
                mpl.patches.Patch(color="red", label="Well"),
                mpl.patches.Patch(color="teal", label="River"),
                mpl.patches.Patch(color="green", label="Drain"),
            ],
            loc="upper left",
        )

        plt.show()
        fig.savefig(paths['figs_path'] / f"{paths['gwf_name']}-head")


def plot_pathpoints(gwf, mf6pl, paths, mp7pl=None, title=None):
    with styles.USGSPlot():
        fig, ax = plt.subplots(ncols=1 if mp7pl is None else 2, nrows=1, figsize=(7, 7))
        if title is not None:
            styles.heading(ax if mp7pl is None else ax[0], heading=title)

        plot_capture_zone(mf6pl)

        plot_points(fig, ax if mp7pl is None else ax[0], gwf, mf6pl)
        if mp7pl is not None:
            plot_points(fig, ax[1], gwf, mp7pl, colorbar=False)

        if mp7pl is not None:
            ax[0].set_xlabel("MODFLOW 6 PRT")
            ax[1].set_xlabel("MODPATH 7")


        plt.show()
        fig.savefig(paths['figs_path'] / f"{paths['sim_name']}-paths-layer.png")


def plot_pathpoints_3d(gwf, mf6pl, paths, sp=1, title=None):
    import pyvista as pv
    from flopy.export.vtk import Vtk

    pv.set_plot_theme("document")
    axes = pv.Axes(show_actor=False, actor_scale=2.0, line_width=5)
    vert_exag = 10
    vtk = Vtk(model=gwf, binary=False, vertical_exageration=vert_exag, smooth=True)
    vtk.add_model(gwf)
    # mf6pl = mf6pl.to_records(index=False)
    vtk.add_pathline_points(mf6pl)
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


    def _plot(paths, screenshot=False):
        p = pv.Plotter(
            window_size=[500, 500],
            off_screen=screenshot,
            notebook=False if screenshot else None,
        )
        p.enable_anti_aliasing()
        if title is not None:
            p.add_title(title, font_size=5)
        p.add_mesh(gwf_mesh, opacity=0.025, style="wireframe")
        p.add_mesh(
            prt_mesh,
            scalars="k" if "k" in prt_mesh.point_data else "ilay",
            # cmap=["red", "red", "green", "blue"],
            point_size=8,
            line_width=3,
            render_points_as_spheres=True,
            render_lines_as_tubes=True,
            smooth_shading=True,
        )

        for key, bc_mesh in bc_meshes.items():

            # plot each of the cells using the standard color
            for bmsh in bc_mesh:
                p.add_mesh(bmsh, color=bc_color_dict[key.upper()], opacity=0.2)

        # p.add_mesh(bed_mesh, color="tan", opacity=0.1)
        # p.remove_scalar_bar()
        p.add_legend(
            labels=[
                ("Well (layer 3)", "red"),
                ("Drain", "green"),
                ("River", "blue"),
            ],
            bcolor="white",
            face="r",
            size=(0.15, 0.15),
        )

        p.camera.zoom(2)
        p.show()
        if screenshot:
            p.screenshot(paths['figs_path'] / f"{paths['sim_name']}-paths-3d.png")

    _plot(paths)
    _plot(paths, screenshot=True)


def plot_all_pathlines(gwf, mf6pl, paths, title=None):
    with styles.USGSPlot():
        fig, ax = plt.subplots(figsize=(7, 7))
        if title is not None:
            styles.heading(ax, heading=title)

        plot_pathlines(ax, gwf, mf6pl, colordest=colordest)

        ax.legend(
            handles=[
                Line2D([0], [0], color="red", lw=1, label="Well"),
                Line2D([0], [0], color="blue", lw=1, label="River"),
                Line2D([0], [0], color="green", lw=1, label="Drain"),
            ],
        )

        plt.show()
        fig.savefig(paths['figs_path'] / f"{paths['sim_name']}-paths.png")


def plot_endpoints(
    gwf, mf6pts, paths, mp7pts=None, title=None, fig_name=None, color="destination"
):
    with styles.USGSPlot():
        fig, ax = plt.subplots(
            ncols=1 if mp7pts is None else 2, nrows=1, figsize=(7, 7)
        )
        if title is not None:
            styles.heading(ax if mp7pts is None else ax[0], heading=title)

        kwargs = {}
        # if color == "destination":
        #     kwargs["colordest"] = colordest

        pts = plot_points(fig, ax if mp7pts is None else ax[0], gwf, mf6pts, **kwargs)
        if mp7pts is not None:
            plot_points(fig, ax[1], gwf, mf6pts, **kwargs)

        if mp7pts is not None:
            ax[0].set_xlabel("MODFLOW 6 PRT")
            ax[1].set_xlabel("MODPATH 7")

        if color == "destination":
            (ax if mp7pts is None else ax[0]).legend(
                handles=[
                    Line2D(
                        [0],
                        [0],
                        color="red",
                        marker="o",
                        markerfacecolor="red",
                        markersize=5,
                        lw=0,
                        label="Well",
                    ),
                    Line2D(
                        [0],
                        [0],
                        color="blue",
                        marker="o",
                        markerfacecolor="blue",
                        markersize=5,
                        lw=0,
                        label="River",
                    ),
                    Line2D(
                        [0],
                        [0],
                        color="green",
                        marker="o",
                        markerfacecolor="green",
                        markersize=5,
                        lw=0,
                        label="Drain",
                    ),
                ],
            )
        else:
            cax = fig.add_axes([0.2, 0.18, 0.6, 0.01])
            cb = plt.colorbar(pts, cax=cax, orientation="horizontal", shrink=0.25)
            cb.set_label("Travel time")

        plt.show()
        fig.savefig(paths['figs_path'] / f"{gwf.name}_plot_endpoints.png")

def plot_capture_zone(pathlines):
    import geopandas as gpd
    from shapely import Polygon
    points = {iprp: [] for iprp in np.unique(pathlines.iprp)}
    for idx, row in pathlines.iterrows():
        points[row.iprp].append((row.x, row.y))
    gdfs = {}
    for key, pnts in points.items():
        p = Polygon(pnts)
        x = p.convex_hull
        gdf = gpd.GeoDataFrame({'name': ['demo'], 'geometry': [x]}, geometry='geometry')
        ax = plt.gca()
        gdf.plot(ax=ax, color='red', alpha=0.2)
        gdfs[key] = gdf


def plot_all(gwf, prt, paths):
    # load results
    head = flopy.utils.HeadFile(os.path.join(gwf.model_ws, gwf.name + ".hds")).get_data()
    mf6pathlines = get_mf6_pathlines(paths['prt_ws'] / paths['trackcsvfile_prt'])
    mp7pathlines = get_mp7_pathlines(
        paths['mp7_ws'] / paths['timeseriesfile_mp7'], paths['mp7_ws'] / paths['endpointfile_mp7'], gwf
    )

    # plot the results
    plot_head(gwf, paths, head=head)
    plot_pathpoints(gwf, mf6pathlines, paths, )  # title="2000-day points, colored by travel time"
    plot_pathpoints_3d(
        gwf, mf6pathlines, paths,  # title="Pathlines, 2000-day points,\ncolored by destination"
    )

    # plot_endpoints(
    #     gwf,
    #     mf6pathlines[(mf6pathlines.ireason == 0) | (mf6pathlines.ireason == 3)],
    #     paths,
    #     title="Release and termination points, colored by destination",
    #     fig_name=f"{paths['sim_name']}-rel-term",
    #     color="destination",
    # )

def get_mf6_pathlines(path):
    # load mf6 pathlines
    pl = pd.read_csv(path)

    # index temporarily by composite key fields
    pl.set_index(["iprp", "irpt", "trelease"], drop=False, inplace=True)

    # # determine which particles ended up in which capture zone
    # pl["destzone"] = pl[pl.istatus > 1].izone
    # pl["dest"] = pl.apply(
    #     lambda row: (
    #         "well"
    #         if (row.destzone == 2 or row.destzone == 3)
    #         else (
    #             "drain"
    #             if row.destzone == 4
    #             else "river"
    #             if row.destzone == 5
    #             else pd.NA
    #         )
    #     ),
    #     axis=1,
    # )
    #
    # # reset index
    # pl.reset_index(drop=True, inplace=True)

    # convert indices to 0-based
    for n in ["imdl", "iprp", "irpt", "ilay", "icell"]:
        pl[n] -= 1

    return pl

def get_mp7_timeseries(path, gwf_model, ref_time_days=0):
    file = flopy.utils.TimeseriesFile(path)
    ts = pd.DataFrame(
        file.get_destination_timeseries_data(list(range(gwf_model.modelgrid.nnodes)))
    )

    # adjust time column since mp7 reports time w.r.t. reference time
    ts["time"] = ts["time"] + ref_time_days

    return ts


def get_mp7_endpoints(path, gwf_model, ref_time_days=0):
    file = flopy.utils.EndpointFile(path)
    ep = pd.DataFrame(
        file.get_destination_endpoint_data(list(range(gwf_model.modelgrid.nnodes)))
    )

    # adjust time column since mp7 reports time w.r.t. reference time
    ep["time"] = ep["time"] + ref_time_days

    return ep

def get_mp7_pathlines(timeseriesfile_path, endpointfile_path, gwf):
    timeseries = get_mp7_timeseries(timeseriesfile_path, gwf, ref_time_days=0)
    endpoints = get_mp7_endpoints(endpointfile_path, gwf, ref_time_days=0)
    return pd.concat([timeseries, endpoints]).reset_index(drop=False)
