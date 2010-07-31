import mpl_toolkits.axisartist.grid_finder as grid_finder
import mpl_toolkits.axisartist.angle_helper as angle_helper

class LocatorDMS(angle_helper.LocatorDMS):
    def set_params(self, **kwargs):
        if "nbins" in kwargs:
            self.den = kwargs.pop("nbins")

class LocatorHMS(angle_helper.LocatorHMS):
    def set_params(self, **kwargs):
        if "nbins" in kwargs:
            self.den = kwargs.pop("nbins")

class FixedLocator(grid_finder.FixedLocator):
    def __init__(self, locs):
        self._locs = locs
        self._factor = None


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



class MaxNLocator(grid_finder.MaxNLocator):
    def set_params(self, **kwargs):
        if "factor" in kwargs:
            self.set_factor(kwargs.pop("factor"))

        grid_finder.MaxNLocator.set_params(self, **kwargs)

