import pyfits

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pywcsgrid2

# read in the image
xray_name="pspc_skyview.fits"
f_xray = pyfits.open(xray_name)

# grid helper to be used.
grid_helper = pywcsgrid2.GridHelper(wcs=f_xray[0].header)
grid_helper.update_wcsgrid_params(label_density=(4,4))

# AxesGrid to display tow images side-by-side
fig = plt.figure(1, (4,3.5))
ax1 = pywcsgrid2.Axes(fig, (0.2, 0.2, 0.75, 0.75), grid_helper=grid_helper)
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
ax1.axis["b=5.5"].label.set_bbox(dict(boxstyle="square,pad=0.1",
                                    fc="w", ec="none", alpha=0.7))
#ax1.axis["b=6"].major_ticklabels.set_visible(False)

ax1.axis["l=82.5"] = ax1["gal"].new_floating_axis(0, 82.5)
ax1["gal"].update_wcsgrid_params(coord_format=("dms","dms"))

plt.draw()
plt.show()

