import pyfits

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.axes_grid import AxesGrid
#from pywcsgrid2.axes_wcs import GridHelperWcs, AxesWcs
import pywcsgrid2

# read in the first image
xray_name="pspc_skyview.fits"
f_xray = pyfits.open(xray_name)
header_xray = f_xray[0].header

# the second image
radio_name="radio_21cm.fits"
f_radio = pyfits.open(radio_name)
header_radio = f_radio[0].header


# grid helper
grid_helper = pywcsgrid2.GridHelper(wcs=header_xray)

# AxesGrid to display tow images side-by-side
fig = plt.figure(1, (6,3.5))
grid = AxesGrid(fig, (0.15, 0.15, 0.8, 0.75), nrows_ncols=(1, 2),
                axes_pad=0.1, share_all=True,
                cbar_mode="each", cbar_location="top", cbar_pad=0,
                axes_class=(pywcsgrid2.Axes, dict(grid_helper=grid_helper)))


ax1 = grid[0]
# use imshow for a simply image display.
im = ax1.imshow(f_xray[0].data, origin="lower", vmin=0., cmap=cm.gray_r,
                interpolation="nearest")
im.set_clim(4.e-05, 0.00018)

ticklocs = [6, 9, 12, 15]

cax1 = grid.cbar_axes[0]
cbar1 = cax1.colorbar(im)
cax1.toggle_label(True)
cax1.set_xticks([t*1.e-5 for t in ticklocs])
cax1.set_xticklabels(["$%d$" % t for t in ticklocs])
#cax1.xaxis.get_major_formatter().set_offset_string(r"$\times 10^{-5}$")
cax1.annotate(r"$\times 10^{-5}$",
              xy=(1,1), xycoords="axes fraction",
              xytext=(0, 15), textcoords="offset points",
              va="bottom", ha="right", size="small")


ax2 = grid[1]
d = f_radio[0].data
# The second image have a different wcs. While imshow works, it will
# interpolate the second image into the image coordinate of the first
# image. You may use pcolormesh when the pixel size of the second
# image is larger than that of the first image. Or you may use
# inshow_affine.

#im2 = ax2[header_radio].pcolormesh(d, cmap=cm.gray_r)
im2 = ax2[header_radio].imshow_affine(d,
                                      cmap=cm.gray_r, origin="lower")

grid.cbar_axes[1].colorbar(im2)
grid.cbar_axes[1].toggle_label(True)

# draw contour. The data points of the contour lines are created in
# the image coordinate of the second image and then are transformed to
# the image coordinate of the first image.
cont = ax2[header_radio].contour(d, colors="k")

# draw contour of the second image in the first axes.
cont2 = ax1[header_radio].contour(d, colors="k")

ax1.add_inner_title("X-ray", loc=2)
ax2.add_inner_title("Radio", loc=2)

ax1.locator_params("both", nbins=2) # since ax1 and ax2 shares a
                                    # grid_helper, it affects not only
                                    # ax1 but also ax2.

plt.show()

