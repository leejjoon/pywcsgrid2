import mpl_toolkits.axisartist.grid_finder as grid_finder
import mpl_toolkits.axisartist.angle_helper as angle_helper
import numpy as np

class LocatorDMS(angle_helper.LocatorDMS):
    def set_params(self, **kwargs):
        if "nbins" in kwargs:
            self.den = int(kwargs.pop("nbins"))

class LocatorHMS(angle_helper.LocatorHMS):
    def set_params(self, **kwargs):
        if "nbins" in kwargs:
            self.den = int(kwargs.pop("nbins"))

class FixedLocator(object):
    def __init__(self, locs, factor=None):
        self._locs = locs
        self._factor = factor


    def __call__(self, v1, v2):
        if self._factor is None:
            v1, v2 = sorted([v1, v2])
        else:
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



class MaxNLocator(grid_finder.MaxNLocator):
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
        
