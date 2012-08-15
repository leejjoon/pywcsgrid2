
import pyfits
import matplotlib.pyplot as plt

fig = plt.figure(1, [5,5])

f = pyfits.open("data/lmc.fits")
d, h = f[0].data, f[0].header

import pywcsgrid2
ax = pywcsgrid2.Axes(fig, [0.2, 0.15, 0.75, 0.8], header=h)
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

