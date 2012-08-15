
import pyfits
import matplotlib.pyplot as plt

fig = plt.figure(1, [5,4.5])

f = pyfits.open("data/lmc.fits")
d, h = f[0].data, f[0].header

import pywcsgrid2
ax = pywcsgrid2.Axes(fig, [0.2, 0.15, 0.65, 0.8], header=h)
fig.add_axes(ax)

im1 = ax.imshow(d, origin="low", vmin=0, vmax=2000,
                cmap=plt.cm.gray_r)

# draw grids
ax.grid()

# grid & ticks in Galactic coordinate.
ax.set_display_coord_system("gal")
ax.set_ticklabel_type("dms", "dms")

# let xaxis display "b", and yaxis "l"
ax.swap_tick_coord()

# turn on top and right ticktables
ax.axis["top","right"].toggle(ticklabels=True)

plt.draw()
plt.show()

