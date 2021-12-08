from pywcsgrid2.axes_wcs import GridHelperWcsFloating, AxesWcs

from mpl_toolkits.axisartist.floating_axes import floatingaxes_class_factory
import matplotlib.axes as maxes

_FloatingAxes = floatingaxes_class_factory(AxesWcs)

class FloatingAxes(_FloatingAxes):
    def __init__(self, fig, rect, header=None, extremes=None, **kwargs):
        if header is None or extremes is None:
            raise ValueError()

        grid_helper = GridHelperWcsFloating(wcs=header, extremes=extremes)
        super(FloatingAxes, self).__init__(fig, rect,
                                           grid_helper=grid_helper,
                                           **kwargs)

        self.set_autoscale_on(False)

    def cla(self):
        super(FloatingAxes, self).cla()

        self.axis["top"].set_ticklabel_direction("+")
        self.axis["top"].set_axislabel_direction("+")
        self.axis["bottom"].set_ticklabel_direction("-")
        self.axis["bottom"].set_axislabel_direction("-")


FloatingSubplot = maxes.subplot_class_factory(FloatingAxes)
