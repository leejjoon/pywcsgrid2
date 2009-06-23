
import pyfits
import matplotlib.pyplot as plt

fig = plt.figure(1, [5,5])

f = pyfits.open("data/lmc.fits")
d, h = f[0].data, f[0].header

from pywcsgrid2.axes_wcs import SubplotWcs
ax = SubplotWcs(fig, 111, header=h)
fig.add_subplot(ax)
fig.subplots_adjust(left=0.2)

im1 = ax.imshow(d, origin="low", vmin=0, vmax=2000,
                cmap=plt.cm.gray_r)

plt.show()

