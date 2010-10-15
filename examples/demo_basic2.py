
import pyfits
import matplotlib.pyplot as plt


f = pyfits.open("data/lmc.fits")
d, h = f[0].data, f[0].header

plt.figure(1, figsize=[5,5])
import pywcsgrid2
ax = pywcsgrid2.axes([0.2, 0.15, 0.7, 0.8], header=h)

im1 = ax.imshow(d, origin="low", vmin=0, vmax=2000,
                cmap=plt.cm.gray_r)

# viewlimits in image coordinate
ax.set_xlim(6, 81)
ax.set_ylim(23, 98)

# draw grids
ax.grid()

# make y ticklabels in absolute degree.
ax.set_ticklabel2_type("absdeg", locs=[-65, -67.5, -70, -72.5])
ax.set_default_label(None, "absdeg") # update the label using the default values

# change grid density of the x axis
ax.locator_params(axis="x", nbins=6)

plt.show()

