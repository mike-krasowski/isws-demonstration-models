from mpl_toolkits.mplot3d import Axes3D
import flopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os

def headContourPlot(model, head, frf, fff,
                    title, subtitle, grid=True, cmap='viridis',
                    plot_disch=True, fs=(8, 8), save_plot=False):
    # data, color-plotting setup
    masked = np.ma.masked_where(model.bas6.ibound.array == 0, head)  # mask inactive cells cells
    ticker = mp.ticker.MaxNLocator(nbins=model.dis.ncol)
    levels = ticker.tick_values(np.min(masked), np.max(masked))  # set contour levels
    cmap = plt.get_cmap(cmap)  # get colormap
    norm = mp.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)  # normalize colormap to dataset levels

    # create figure and plots
    fig, ax = plt.subplots(figsize=fs)  # create figure
    modelmap = flopy.plot.PlotMapView(model=model, layer=0, ax=ax)  # use modelmap to attach plot to model
    pcolor = modelmap.plot_array(masked, norm=norm, animated=True)  # create contours
    if grid is True:
        grid = modelmap.plot_grid()  # plot grid if True
    if plot_disch is True:
        modelmap.plot_vector(frf[0], fff[0])  # create discharge arrows if true

    # figure display parameters
    plt.xlabel('Lx (m)', fontsize=13)
    plt.ylabel('Ly (m)', fontsize=13)
    plt.suptitle(title, fontsize=20, fontweight='bold', x=.44, y=.98)  # bold title
    plt.title(subtitle, fontsize=14, y=1.015)  # subtitle
    cb = plt.colorbar(pcolor)
    cb.set_label('Head (m)', fontsize=13, labelpad=12)  # add colorbar for head data
    plt.show()
    return (fig)


def headSurfacePlot(model, Lx, Ly, head,
                    title='', cmap='viridis', fs=(12, 5)):
    # create 3d figure
    fig_3d = plt.figure(figsize=fs)
    ax = fig_3d.add_subplot(111, projection='3d')

    # set up data space
    #    masked = np.where(model.bas6.ibound.array !=0, head, np.nan)
    masked = head[:, 1:-1, :]
    x = np.linspace(0, Lx, model.dis.ncol)  # set x domain using model dis pkg
    y = np.linspace(0, Ly, model.dis.nrow)  # set y domain using model dis pkg
    #    x, y = np.meshgrid(x, y) #mesh domain
    x, y = np.meshgrid(x, y[1:-1])  # mesh domain

    # set color levels
    levels = np.linspace(np.nanmin(masked),
                         np.nanmax(masked),
                         model.dis.ncol - 1)
    cmap = plt.get_cmap(cmap)  # get colormap
    norm = mp.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)  # normalize colormap to dataset values
    surf = ax.plot_surface(x, y, np.flipud(masked[0]),
                           cmap=cmap, linewidth=0, antialiased=False,
                           label='head', norm=norm)  # plot active head surface
    surf.set_facecolor((0, 0, 0, 0))
    cb = fig_3d.colorbar(surf, shrink=0.5, aspect=5)
    cb.set_label('Head (m)', fontsize=15, labelpad=10)  # add colorbar
    ax.set_xlabel('Lx (m)', fontsize=13, labelpad=10)
    ax.set_ylabel('Ly (m)', fontsize=13, labelpad=10)
    ax.set_title(title, fontsize=15, y=1.05)
    # suppress matplotlib warning (3d plot is angry at nan values but still works)
    import warnings
    warnings.filterwarnings("ignore")
    plt.show()
    return (fig_3d)