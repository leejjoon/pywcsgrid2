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


# grid helper to be used.
grid_helper = pywcsgrid2.GridHelperSky(wcs=header_xray)

# AxesGrid to display tow images side-by-side
fig = plt.figure(1, (6,3.5))
grid = AxesGrid(fig, (0.15, 0.15, 0.8, 0.75), nrows_ncols=(1, 2),
                axes_pad=0.1, share_all=True,
                axes_class=(pywcsgrid2.Axes, dict(grid_helper=grid_helper)))


ax1 = grid[0]

#from mpl_toolkits.axes_grid1.angle_helper import LocatorDMS
#from pywcsgrid2.axes_wcs import FormatterDMSDelta

#ax1.get_grid_helper().set_ticklabel_mode("delta", center_pixel=(50, 50))
ax1.set_ticklabel_type("delta",
                       center_pixel=(146, 139))
ax1.locator_params(nbins=3)


# use imshow for a simply image display.
im = ax1.imshow(f_xray[0].data, origin="lower", vmin=0., cmap=cm.gray_r,
                interpolation="nearest")
im.set_clim(4.e-05, 0.00018)


ax2 = grid[1]
d = f_radio[0].data
# The second image have a different wcs. While imshow works (it
# will interpolate the second image into the image coordinate of the
# first image), pcolormesh is prefered when the pixel size of the
# second image is larger than that of the first image.
im2 = ax2[header_radio].pcolormesh(d, cmap=cm.gray_r)

# draw contour. The data points of the contour lines are created in
# the image coordinate of the second image and then are transformed to
# the image coordinate of the first image.
cont = ax2[header_radio].contour(d, colors="k")

# draw contour of the second image in the first axes.
cont2 = ax1[header_radio].contour(d, colors="k")

ax1.set_title("X-ray")
ax2.set_title("Radio")

plt.draw()
plt.show()

