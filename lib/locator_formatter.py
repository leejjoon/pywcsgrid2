from matplotlib import ticker as mticker
import mpl_toolkits.axisartist.grid_finder as grid_finder
import mpl_toolkits.axisartist.angle_helper as angle_helper
import numpy as np

from mpl_toolkits.axisartist.angle_helper import LocatorH, LocatorHM, LocatorHMS, \
     LocatorD, LocatorDM, LocatorDMS

class FixedLocator(object):
    def __init__(self, locs, factor=1):
        self._locs = locs
        self._factor = factor


    def __call__(self, v1, v2):
        v1, v2 = sorted([v1*self._factor, v2*self._factor])
        locs = np.array([l for l in self._locs if ((v1 <= l) and (l <= v2))])

        return locs, len(locs), self._factor

    def set_factor(self, f):
        self._factor = f

    def set_params(self, **kwargs):
        if "factor" in kwargs:
            self.set_factor(kwargs.pop("factor"))
        if 'locs' in kwargs:
            self._locs = kwargs.pop('locs')


# This was originally part of the AxisArtist. But set_factor method was removed
# in mpl3.5. We simply copy the orinal MaxNLocator here.
class MaxNLocator(mticker.MaxNLocator):
    def __init__(self, nbins=10, steps=None,
                 trim=True,
                 integer=False,
                 symmetric=False,
                 prune=None):
        # trim argument has no effect. It has been left for API compatibility
        super().__init__(nbins, steps=steps, integer=integer,
                         symmetric=symmetric, prune=prune)
        self.create_dummy_axis()
        self._factor = 1

    def __call__(self, v1, v2):
        locs = super().tick_values(v1 * self._factor, v2 * self._factor)
        return np.array(locs), len(locs), self._factor

    def set_factor(self, f):
        self._factor = f

    def set_params(self, **kwargs):
        if "factor" in kwargs:
            self.set_factor(kwargs.pop("factor"))

        grid_finder.MaxNLocator.set_params(self, **kwargs)


class FixedFormatter(object):
    def __init__(self, labels):
        """
        format_dict : dictionary for format strings to be used.
        formatter : fall-back formatter
        """
        super(FixedFormatter, self).__init__()
        self._labels = labels

    def __call__(self, direction, factor, values):
        """
        factor is ignored.
        """

        return self._labels

import mpl_toolkits.axisartist.grid_finder as grid_finder
class FormatterPrettyPrint(grid_finder.FormatterPrettyPrint):
    def __init__(self, useMathText=True):
        grid_finder.FormatterPrettyPrint.__init__(self)
        #self._fmt._useMathText = useMathText
        self._fmt._usetex = True


