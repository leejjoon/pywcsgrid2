import pyfits
import pywcs

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid.axes_grid import AxesGrid
from pywcsgrid2.axes_wcs import GridHelperWcs, AxesWcs

# read in the image
xray_name="pspc_skyview.fits"
f_xray = pyfits.open(xray_name)
wcs_xray = pywcs.WCS(f_xray[0].header)

# grid helper to be used.
grid_helper = GridHelperWcs(wcs=wcs_xray)
grid_helper.update_wcsgrid_params(label_density=(2,2))

# AxesGrid to display tow images side-by-side
fig = plt.figure(1, (4,3.5))
ax1 = AxesWcs(fig, (0.2, 0.2, 0.75, 0.75), grid_helper=grid_helper)
fig.add_axes(ax1)

# use imshow for a simply image display.
im = ax1.imshow(f_xray[0].data, origin="lower", vmin=0., cmap=cm.gray_r,
                interpolation="nearest")
im.set_clim(4.e-05, 0.00018)


ax1.axis["b=5"] = axis = ax1["gal"].get_grid_helper().new_floating_axis(1, 5, axes=ax1)
axis.label.set_text(r"$b=5^{\circ}$")
axis.major_ticklabels.set_visible(False)

def grow(im, dpi):
    "enlarge the area"
    import numpy as np
    import scipy.ndimage as NI

    pad = 5
    ny, nx, depth = im.shape
    new_alpha = np.zeros([pad*2+ny, pad*2+nx], dtype="d")
    new_alpha[pad:-pad, pad:-pad] = im[:,:,-1]
    #new_alpha[pad:-pad, pad:-pad] = (1.-im[:,:,0])
    alpha2 = NI.grey_dilation(new_alpha, size=(3, 3))
    new_im = np.zeros([pad*2+ny, pad*2+nx, depth], dtype="d")
    new_im[:,:,-1] = alpha2
    new_im[:,:,:-1] = 1.
    offsetx, offsety = -pad, -pad
    return new_im, offsetx, offsety

from matplotlib.artist import Artist

class FilteredArtistList(Artist):
    """
    A simple container to draw filtered artist.
    """
    def __init__(self, artist_list, filter):
        self._artist_list = artist_list
        self._filter = filter
        Artist.__init__(self)
        
    def draw(self, renderer):
        renderer.start_rasterizing()
        renderer.start_filter()
        for a in self._artist_list:
            a.draw(renderer)
        renderer.stop_filter(self._filter)
        renderer.stop_rasterizing()



white_glows = FilteredArtistList([axis], grow)
ax1.add_artist(white_glows)
white_glows.set_zorder(axis.get_zorder()-0.1)

plt.draw()
plt.show()

