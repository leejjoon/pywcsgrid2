
import numpy as np


from matplotlib.patches import FancyArrowPatch
from matplotlib.text import Text

from matplotlib.offsetbox import AnchoredOffsetbox, AuxTransformBox

from pywcsgrid2.wcs_helper import estimate_angle


class AnchoredCompass(AnchoredOffsetbox):
    def __init__(self, ax, transSky2Pix, loc,
                 arrow_fraction=0.15,
                 txt1="E", txt2="N",
                 delta_a1=0, delta_a2=0,
                 pad=0.1, borderpad=0.5, prop=None, frameon=False,
                 ):
        """
        Draw an arrows pointing the directions of E & N

        arrow_fraction : length of the arrow as a fraction of axes size

        pad, borderpad in fraction of the legend font size (or prop)
        """

        self._ax = ax
        self._transSky2Pix = transSky2Pix
        self._box = AuxTransformBox(ax.transData)
        self.delta_a1, self.delta_a2 = delta_a1, delta_a2
        self.arrow_fraction = arrow_fraction

        kwargs = dict(mutation_scale=11,
                      shrinkA=0,
                      shrinkB=5)

        self.arrow1 = FancyArrowPatch(posA=(0, 0), posB=(1, 1),
                                      arrowstyle="->",
                                      arrow_transmuter=None,
                                      connectionstyle="arc3",
                                      connector=None,
                                      **kwargs)
        self.arrow2 = FancyArrowPatch(posA=(0, 0), posB=(1, 1),
                                      arrowstyle="->",
                                      arrow_transmuter=None,
                                      connectionstyle="arc3",
                                      connector=None,
                                      **kwargs)


        x1t, y1t, x2t, y2t = 1, 1, 1, 1
        self.txt1 = Text(x1t, y1t, txt1, rotation=0,
                         rotation_mode="anchor",
                         va="center", ha="right")
        self.txt2 = Text(x2t, y2t, txt2, rotation=0,
                         rotation_mode="anchor",
                         va="bottom", ha="center")


        self._box.add_artist(self.arrow1)
        self._box.add_artist(self.arrow2)

        self._box.add_artist(self.txt1)
        self._box.add_artist(self.txt2)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=self._box,
                                   prop=prop,
                                   frameon=frameon)

    def set_path_effects(self, path_effects):
        for a in [self.arrow1, self.arrow2, self.txt1, self.txt2]:
            a.set_path_effects(path_effects)

    def _update_arrow(self, renderer):
        ax = self._ax

        x0, y0 = ax.viewLim.x0, ax.viewLim.y0
        a1, a2 = estimate_angle(self._transSky2Pix, x0, y0)
        a1, a2 = a1+self.delta_a1, a2+self.delta_a2

        D = min(ax.viewLim.width, ax.viewLim.height)
        d = D * self.arrow_fraction
        x1, y1 = x0+d*np.cos(a1/180.*np.pi), y0+d*np.sin(a1/180.*np.pi)
        x2, y2 = x0+d*np.cos(a2/180.*np.pi), y0+d*np.sin(a2/180.*np.pi)

        self.arrow1.set_positions((x0, y0), (x1, y1))
        self.arrow2.set_positions((x0, y0), (x2, y2))

        d2 = d
        x1t, y1t = x0+d2*np.cos(a1/180.*np.pi), y0+d2*np.sin(a1/180.*np.pi)
        x2t, y2t = x0+d2*np.cos(a2/180.*np.pi), y0+d2*np.sin(a2/180.*np.pi)

        self.txt1.set_position((x1t, y1t))
        self.txt1.set_rotation(a1-180)
        self.txt2.set_position((x2t, y2t))
        self.txt2.set_rotation(a2-90)


    def draw(self, renderer):
        self._update_arrow(renderer)
        super(AnchoredCompass, self).draw(renderer)
