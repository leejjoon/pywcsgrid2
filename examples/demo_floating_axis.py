import pyfits

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pywcsgrid2
import matplotlib.patheffects as patheffects


# read in the image
xray_name="pspc_skyview.fits"
f_xray = pyfits.open(xray_name)

fig = plt.figure(1, (4,3.5))
ax1 = pywcsgrid2.Axes(fig, (0.2, 0.2, 0.75, 0.75),
                      header=f_xray[0].header)
fig.add_axes(ax1)

# use imshow for a simple image display.
im = ax1.imshow(f_xray[0].data, origin="lower", vmin=0.,
                cmap=cm.gray_r,
                interpolation="nearest")
im.set_clim(4.e-05, 0.00018)


# create a new floating axis in Galactic coordinate
#grid_helper_gal = ax1["gal"].get_grid_helper()
ax1.axis["b=5.5"] = ax1["gal"].new_floating_axis(1, 5.5)
ax1.axis["b=5.5"].label.set_text(r"$b=5.5^{\circ}$")

# use path_effects
my_path_effects = [patheffects.withStroke(linewidth=3,
                                          foreground="w")]
ax1.axis["b=5.5"].major_ticklabels.set_path_effects(my_path_effects)
ax1.axis["b=5.5"].label.set_path_effects(my_path_effects)

ax1.axis["l=82.5"] = ax1["gal"].new_floating_axis(0, 82.5)

plt.draw()
plt.show()


