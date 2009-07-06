import matplotlib.pyplot as plt
import pyfits
from mpl_toolkits.axes_grid.axes_grid import AxesGrid
import pywcs
from pywcsgrid2.axes_wcs import AxesWcs, GridHelperWcs
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes, mark_inset
from mpl_toolkits.axes_grid.axislines import Axes

if 1:

    f = pyfits.open("data/E0102_E0300-4000.fits")
    d = f[0].data
    h = f[0].header

    fig = plt.figure(1)
    fig.clf()

    wcs = pywcs.WCS(h)
    grid_helper = GridHelperWcs(wcs)
    grid_helper.update_wcsgrid_params(label_density=(5,5))

    grid = AxesGrid(fig, (1, 1, 1), nrows_ncols = (1, 1),
                    cbar_mode="single", cbar_pad="2%",
                    cbar_location="right",
                    axes_class=(AxesWcs, dict(grid_helper=grid_helper)))


    im = grid[0].imshow(d, origin="lower", cmap=plt.cm.gray_r,
                        interpolation="nearest")

    grid.cbar_axes[0].colorbar(im)
    grid.cbar_axes[0].axis["right"].major_ticklabels.set_visible(True)

    #grid_helper_inset = GridHelperWcs(wcs)
    #grid_helper_inset.update_wcsgrid_params(label_density=(3,3))

    axins = zoomed_inset_axes(grid[0], 3, loc=1,
                              axes_class=Axes)
    #                          axes_kwargs=dict(grid_helper=grid_helper_inset))

    for a in axins.axis.values():
        a.major_ticklabels.set_visible(False)
        a.major_ticks.set_visible(False)

    axins.imshow(d, origin="lower", interpolation="nearest",
                 cmap=plt.cm.gray_r)

    axins.set_xlim(23, 30)
    axins.set_ylim(22, 29)

    mark_inset(grid[0], axins, loc1=2, loc2=4, fc="none", ec="0.5")

    plt.draw()
    plt.show()

