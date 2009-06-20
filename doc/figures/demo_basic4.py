
import pyfits
import pywcs
import matplotlib.pyplot as plt

fig = plt.figure(1, [5,5])

f = pyfits.open("data/lmc.fits")
d, h = f[0].data, f[0].header
wcs = pywcs.WCS(h)

from pywcsgrid2.axes_wcs import AxesWcs
ax = AxesWcs(fig, [0.2, 0.15, 0.75, 0.8], wcs=wcs)
fig.add_axes(ax)

im1 = ax.imshow(d, origin="low", vmin=0, vmax=2000,
                cmap=plt.cm.gray_r)

ax.set_display_coord_system("fk5")

# draw grids
ax.grid()

# (alpha, delta) in degree
ax["fk4"].plot([x/24.*360 for x in [4, 5, 6]],
               [-74, -70, -66], "ro-")

# (l, b)  in degree
ax["gal"].plot([(285), (276.)],
               [(-30), (-36)])


plt.draw()
plt.show()

