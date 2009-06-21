
import pyfits
import pywcs
import matplotlib.pyplot as plt

fig = plt.figure(1, [5,5])

f = pyfits.open("data/lmc.fits")
d, h = f[0].data, f[0].header
wcs = pywcs.WCS(h)

from pywcsgrid2.axes_wcs import AxesWcs
ax = AxesWcs(fig, [0.2, 0.15, 0.7, 0.8], wcs=wcs)
fig.add_axes(ax)

im1 = ax.imshow(d, origin="low", vmin=0, vmax=2000,
                cmap=plt.cm.gray_r)

# viewlimits in image coordinate
ax.set_xlim(6, 81)
ax.set_ylim(23, 98)

#ax.axis["bottom"].label.set_text(r"$\alpha_{1950}$")
#ax.axis["left"].label.set_text(r"$\delta_{1950}$")

# draw grids
ax.grid()

# change grid density
ax.get_grid_helper().update_wcsgrid_params(label_density=[6,4])

plt.draw()
plt.show()

