
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

plt.show()

