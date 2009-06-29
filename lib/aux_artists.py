
import numpy as np


from matplotlib.patches import FancyArrowPatch
from matplotlib.text import Text

from matplotlib.offsetbox import AnchoredOffsetbox, AuxTransformBox

from pywcsgrid2.wcs_helper import estimate_angle


class AnchoredCompass(AnchoredOffsetbox):
    def __init__(self, ax, transSky2Pix, loc,
                 arrow_length=0.15,
                 txt1="E", txt2="N",
                 delta_a1=0, delta_a2=0,
                 pad=0.1, borderpad=0.5, prop=None, frameon=False,
                 ):
        """
        Draw an arrows pointing the directions of E & N

        arrow_length : length of the arrow as a fraction of axes size

        pad, borderpad in fraction of the legend font size (or prop)
        """
        self._box = AuxTransformBox(ax.transData)
        #self.ellipse = Ellipse((0,0), width, height, angle)


        #wcs_proj = axes.projection
        x0, y0 = ax.viewLim.x0, ax.viewLim.y0
        a1, a2 = estimate_angle(transSky2Pix, x0, y0)
        a1, a2 = a1+delta_a1, a2+delta_a2

        D = min(ax.viewLim.width, ax.viewLim.height)
        d = D * arrow_length
        x1, y1 = x0+d*np.cos(a1/180.*np.pi), y0+d*np.sin(a1/180.*np.pi)
        x2, y2 = x0+d*np.cos(a2/180.*np.pi), y0+d*np.sin(a2/180.*np.pi)

        kwargs = dict(mutation_scale=11,
                      shrinkA=0,
                      shrinkB=5)
        arrow_E = FancyArrowPatch(posA=(x0, y0), posB=(x1, y1),
                                  arrowstyle="->",
                                  arrow_transmuter=None,
                                  connectionstyle="arc3",
                                  connector=None,
                                  **kwargs)
        arrow_N = FancyArrowPatch(posA=(x0, y0), posB=(x2, y2),
                                  arrowstyle="->",
                                  arrow_transmuter=None,
                                  connectionstyle="arc3",
                                  connector=None,
                                  **kwargs)

        d2 = d
        x1t, y1t = x0+d2*np.cos(a1/180.*np.pi), y0+d2*np.sin(a1/180.*np.pi)
        x2t, y2t = x0+d2*np.cos(a2/180.*np.pi), y0+d2*np.sin(a2/180.*np.pi)


        txt1 = Text(x1t, y1t, txt1, rotation=a1-180,
                    rotation_mode="anchor",
                    va="center", ha="right")
        txt2 = Text(x2t, y2t, txt2, rotation=a2-90,
                    rotation_mode="anchor",
                    va="bottom", ha="center")


        self._box.add_artist(arrow_E)
        self._box.add_artist(arrow_N)

        self._box.add_artist(txt1)
        self._box.add_artist(txt2)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=self._box,
                                   prop=prop,
                                   frameon=frameon)

