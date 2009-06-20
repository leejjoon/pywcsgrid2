
import pyfits
import pywcs
import matplotlib.pyplot as plt

fig = plt.figure(1, [5,4.5])

f = pyfits.open("data/lmc.fits")
d, h = f[0].data, f[0].header
wcs = pywcs.WCS(h)

from pywcsgrid2.axes_wcs import AxesWcs
ax = AxesWcs(fig, [0.2, 0.15, 0.65, 0.8], wcs=wcs)
fig.add_axes(ax)

im1 = ax.imshow(d, origin="low", vmin=0, vmax=2000,
                cmap=plt.cm.gray_r)

# draw grids
ax.grid()

# grid & ticks in Galactic coordinate.
ax.set_display_coord_system("gal")

# let xaxis display "b", and yaxis "l"
ax.swap_tick_coord()

# turn on top and right ticktables
ax.axis["top"].major_ticklabels.set_visible(True)
ax.axis["right"].major_ticklabels.set_visible(True)

plt.draw()
plt.show()

