
import matplotlib.axes as maxes
from matplotlib.cbook import is_string_like

import numpy as np
from pywcsgrid2.wcs_transforms import WcsSky2PixelTransform, WcsPixel2SkyTransform,\
     WcsSky2SkyTransform

from matplotlib import rcParams

import matplotlib.collections as mcoll
def rasterized_draw(self, renderer, *kl, **kwargs):
    renderer.start_rasterizing()
    self.draw_orig(renderer, *kl, **kwargs)
    renderer.stop_rasterizing()
mcoll.QuadMesh.draw_orig = mcoll.QuadMesh.draw
mcoll.QuadMesh.draw = rasterized_draw


#from  axislines import AxisLineHelper, GridHelperRectlinear, AxisArtist
from  mpl_toolkits.axes_grid.axislines import GridlinesCollection

#import pywcs
from kapteyn_helper import get_kapteyn_projection
from pywcsgrid2.kapteyn_helper import coord_system_guess

from mpl_toolkits.axes_grid.parasite_axes import HostAxes, ParasiteAxesAuxTrans
import mpl_toolkits.axes_grid.grid_finder as grid_finder

GridFinderBase = grid_finder.GridFinderBase



from mpl_toolkits.axes_grid.angle_helper import ExtremeFinderCycle
from mpl_toolkits.axes_grid.angle_helper import LocatorDMS, LocatorHMS, FormatterDMS, FormatterHMS

from mpl_toolkits.axes_grid.grid_helper_curvelinear \
     import GridHelperCurveLinear, FixedAxisArtistHelper

import weakref
class WcsTransCollection(object):
    def __init__(self):
        self._wcs_trans_collection = dict()

    def _get_wcs_hash(self, wcs):
        return wcs

    def _get_wcs_trans_dict(self, wcs):
        wcs_hash = self._get_wcs_hash(wcs)
        if wcs_hash in self._wcs_trans_collection:
            return self._wcs_trans_collection[wcs_hash]
        else:
            wcs_trans_dict = weakref.WeakValueDictionary()
            wcs_trans_dict[None] = WcsSky2PixelTransform(wcs)
            self._wcs_trans_collection[wcs_hash] = wcs_trans_dict

        return wcs_trans_dict

    def get_wcs_trans(self, wcs, coord_system):
        wcs_trans_dict = self._get_wcs_trans_dict(wcs)
        if coord_system in wcs_trans_dict:
            return wcs_trans_dict[coord_system]
        else:
            wcs_trans = WcsSky2PixelTransform(wcs,
                                              src_coord=coord_system)
            wcs_trans_dict[coord_system] = wcs_trans
            return wcs_trans

wcs_trans_collection = WcsTransCollection()
# parasite axes needs to be modified to use wcs_trans_collection

class GridHelperWcs(GridHelperCurveLinear):

    WCS_TRANS_COLLECTION = wcs_trans_collection

    @classmethod
    def get_wcs_trans(cls, wcs, coord):
        return cls.WCS_TRANS_COLLECTION.get_wcs_trans(wcs, coord)

    def __init__(self, wcs, orig_coord=None):

        #self.header = header
        #wcs = pywcs.WCS(header)
        projection = get_kapteyn_projection(wcs)
        self.projection = projection

        if orig_coord is None:
            ctype1, ctype2 = self.projection.ctype
            equinox = self.projection.equinox
            coord_guess = coord_system_guess(ctype1, ctype2, equinox)

            #coord_guess = coord_system_guess(wcs.wcs.ctype[0],
            #                                 wcs.wcs.ctype[1],
            #                                 wcs.wcs.equinox)
            if coord_guess is None:
                raise ValueError("Unknown coordinate system, %, %s, equinox=%.1f" \
                                 % (ctype1, ctype2,
                                    equinox))
            else:
                self._wcsgrid_orig_coord_system = coord_guess
        else:
            self._wcsgrid_orig_coord_system = orig_coord

        self.wcsgrid = None
        self._old_values = None

        _wcs_trans = self.get_wcs_trans(self.projection, None)
        self._wcsgrid_display_coord_system = None

        self._wcsgrid_params = dict(coord_format=("hms", "dms"),
                                    label_density=(4, 4),
                                    grid1=[], grid2=[])

        GridHelperCurveLinear.__init__(self, _wcs_trans,
                                       extreme_finder=ExtremeFinderCycle(20,20),
                                       grid_locator1=LocatorHMS(4),
                                       grid_locator2=LocatorDMS(4),
                                       tick_formatter1=FormatterHMS(),
                                       tick_formatter2=FormatterDMS())


    def tr(self, lon, lat):
        return self.projection.topixel((lon, lat))
        #ll = np.concatenate((lon[:,np.newaxis], lat[:,np.newaxis]), 1)
        #origin=0
        #xy = self.wcs.wcs.s2p(ll, origin)["pixcrd"]
        #return xy[:,0], xy[:,1]

    def inv_tr(self, x, y):
        return self.projection.toworld((x, y))
        #xy = np.concatenate((x[:,np.newaxis], y[:,np.newaxis]),
        #                    1)
        #origin=0
        #ll = self.wcs.wcs.p2s(xy, origin)['world']
        #return ll[:,0], ll[:,1]


    def update_wcsgrid_params(self, **kwargs):
        """
        coord_format="GE",
        label_density=(6, 6),
        grid1=[], grid2=[]):
        """

        self._wcsgrid_params.update(**kwargs)

        f1, f2 = self._wcsgrid_params["coord_format"]
        locator_selector = dict(hms=LocatorHMS,
                                dms=LocatorDMS).__getitem__
        formatter_selector = dict(hms=FormatterHMS,
                                  dms=FormatterDMS).__getitem__

        locator1 = locator_selector(f1)
        locator2 = locator_selector(f2)
        formatter1 = formatter_selector(f1)
        formatter2 = formatter_selector(f2)

        lx, ly = self._wcsgrid_params["label_density"]
        grid_locator1=locator1(lx)
        grid_locator2=locator2(ly)
        self.update_grid_finder(grid_locator1=grid_locator1,
                                grid_locator2=grid_locator2,
                                tick_formatter1=formatter1(),
                                tick_formatter2=formatter2(),
                                )



    def get_wcsgrid_params(self):
        return self._wcsgrid_params


    def set_display_coord_system(self, coord_system):
        if self._wcsgrid_orig_coord_system is None:
            raise ValueError("_pywcsgrid_orig_coord_system is not set")

        self._wcsgrid_display_coord_system = coord_system
        _wcs_trans = self.get_wcs_trans(self.projection, coord_system)

        self.update_grid_finder(aux_trans=_wcs_trans)


    def get_display_coord_system(self):
        return self._wcsgrid_display_coord_system



class ParasiteAxesSky(ParasiteAxesAuxTrans):
    def __init__(self, parent, src_coord, *kl, **kwargs):

        dest_coord = parent.get_grid_helper()._wcsgrid_orig_coord_system
        tr = WcsSky2SkyTransform(src_coord, dest_coord) + \
             parent._wcsgrid_wcsaxes[0].transAux


        grid_helper = GridHelperWcs(parent.projection)
        grid_helper.set_display_coord_system(src_coord)

        kwargs["grid_helper"] = grid_helper

        ParasiteAxesAuxTrans.__init__(self, parent, tr, *kl, **kwargs)


        self.src_coord = src_coord
        self.dest_coord = dest_coord


    def _get_gridlines(self):
        gridlines = GridlinesCollection(None,
                                        transform=self._parent_axes.transData,
                                        colors=rcParams['grid.color'],
                                        linestyles=rcParams['grid.linestyle'],
                                        linewidths=rcParams['grid.linewidth'])
        self._parent_axes._set_artist_props(gridlines)
        gridlines.set_grid_helper(self.get_grid_helper())

        return gridlines

    def _init_axislines(self):
        self._axislines = self.AxisDict(self)
        new_fixed_axis = self.get_grid_helper().new_fixed_axis
        for loc in ["bottom", "top", "left", "right"]:
            self._axislines[loc] = new_fixed_axis(loc=loc,
                                                  axes=self._parent_axes)

        for axisline in [self._axislines["top"], self._axislines["right"]]:
            axisline.label.set_visible(False)
            axisline.major_ticklabels.set_visible(False)
            axisline.minor_ticklabels.set_visible(False)


    def new_floating_axis(self, nth_coord, value,
                          tick_direction="in",
                          label_direction="top",
                          ):
        gh = self.get_grid_helper()
        axis = gh.new_floating_axis(nth_coord, value,
                                    tick_direction=tick_direction,
                                    label_direction=label_direction,
                                    axes=self._parent_axes)
        return axis


    def update_wcsgrid_params(self, **ka):
        self.get_grid_helper().update_wcsgrid_params(**ka)
        


def get_transformed_image(Z, tr, extent=None, oversample=1.5):

    from scipy.ndimage import map_coordinates

    ny, nx = Z.shape
    if extent is None:
        #extent = 0, nx-1, 0, ny-1
        x1, x2, y1, y2 = 0, nx-1, 0, ny-1
    else:
        x1, x2, y1, y2 = extent

    transformed_edge = tr.transform([(x1, y1), (x1, y2), (x2, y1), (x2, y2)])
    X1 = transformed_edge[:,0].min()
    X2 = transformed_edge[:,0].max()
    Y1 = transformed_edge[:,1].min()
    Y2 = transformed_edge[:,1].max()

    nX, nY = int(nx*oversample), int(ny*oversample),
    mX, mY = np.meshgrid(np.linspace(X1, X2, nX),
                         np.linspace(Y1, Y2, nY))

    mXY = np.concatenate([mX[:,:,np.newaxis],
                          mY[:,:,np.newaxis]],2).reshape(nX*nY, 2)

    tXY = tr.inverted().transform(mXY)

    O = map_coordinates(Z, [tXY[:,1].reshape(nY, nX),
                            tXY[:,0].reshape(nY, nX)],
                        mode='constant', cval=np.nan,
                        prefilter=True, order=0)

    o_extent = X1, X2, Y1, Y2
    return O, o_extent


class ParasiteAxesWcs(ParasiteAxesAuxTrans):
    def __init__(self, parent, wcs, *kl, **kwargs):

        dest_coord = parent.get_grid_helper()._wcsgrid_orig_coord_system
        tr = WcsPixel2SkyTransform(wcs, dest_coord) + \
             parent._wcsgrid_wcsaxes[0].transAux
        ParasiteAxesAuxTrans.__init__(self, parent, tr, *kl, **kwargs)

    def imshow(self, Z, oversample=1.5, sub_range=None, **kwargs):
        try:
            from scipy.ndimage import map_coordinates
        except ImportError:
            raise ImportError("this routine requires scipy installed. You may use pcolormesh instead but will be much slower.")

        if "extent" in kwargs:
            raise ValueError("extent is not allowed")

        O, oe  = get_transformed_image(Z, self.transAux,
                                       extent=sub_range, oversample=oversample)
        Om = np.ma.array(O, mask=np.isnan(O))

        ax = self._parent_axes
        aox = ax.get_autoscalex_on()
        aoy = ax.get_autoscaley_on()
        ax.set_autoscalex_on(False)
        ax.set_autoscaley_on(False)

        im = ax.imshow(Om, extent=oe, **kwargs)

        ax.set_autoscalex_on(aox)
        ax.set_autoscaley_on(aoy)

        return im



class AxesWcs(HostAxes):

    def __init__(self, *kl, **kw):

        if "grid_helper" not in kw:
            header = kw.pop("header", None)
            wcs = kw.pop("wcs", None)
            if (header is not None) and (wcs is None):
                #self._init_kapteyn_projection(header)
                #self._wcs = pywcs.WCS(header)
                self.projection = get_kapteyn_projection(header)
            elif (header is None) and (wcs is not None):
                #self._init_kapteyn_projection(wcs)
                #self._wcs = wcs
                self.projection = get_kapteyn_projection(wcs)
            else:
                raise ValueError("wcs")

            #grid_helper = GridHelperWcs(self._wcs)
            grid_helper = GridHelperWcs(self.projection)
            kw["grid_helper"] = grid_helper
        else:
            self.projection = kw["grid_helper"].projection

        super(AxesWcs, self).__init__(*kl, **kw)

        self._init_parasites()

        self.set_default_label()


    def _init_parasites(self):
        ax = ParasiteAxesAuxTrans(self,
                                  #WcsSky2PixelTransform(self._wcs),
                                  WcsSky2PixelTransform(self.projection),
                                  viewlim_mode="equal")
        self._wcsgrid_wcsaxes = {0:ax}
        self.parasites.append(ax)


    def __getitem__(self, key):

        # check if key is a valid coord_sys instance
        if key == 0: # or isinstance(key, _coord_sys_dict) :
            pass
        elif is_string_like(key):
            pass
            #key = _coord_sys_dict[key.lower()]
        #elif isinstance(key, pywcs.WCS):
        #    pass
        else:
            try:
                key = get_kapteyn_projection(key)
            except:
                raise ValueError("invalide key : %s" % repr(key))

        if key not in self._wcsgrid_wcsaxes:
            if self.get_grid_helper()._wcsgrid_orig_coord_system == key:
                self._wcsgrid_wcsaxes[key] = self._wcsgrid_wcsaxes[0]
            else:
                orig_coord = self.get_grid_helper()._wcsgrid_orig_coord_system

                if is_string_like(key):
                    ax = ParasiteAxesSky(self, key)
                #if isinstance(key, pywcs.WCS):
                #    ax = ParasiteAxesWcs(self, key)
                else:
                    ax = ParasiteAxesWcs(self, key)
                    #ax = ParasiteAxesSky(self, key)

                self._wcsgrid_wcsaxes[key] = ax
                self.parasites.append(ax)


        return self._wcsgrid_wcsaxes[key]


    def update_wcsgrid_params(self, **ka):
        self.get_grid_helper().update_wcsgrid_params(**ka)
        

    def set_display_coord_system(self, c):
        self.get_grid_helper().set_display_coord_system(c)

        self.axis["bottom"].get_helper().change_tick_coord(0)
        self.axis["top"].get_helper().change_tick_coord(0)
        self.axis["left"].get_helper().change_tick_coord(1)
        self.axis["right"].get_helper().change_tick_coord(1)

        self.set_default_label()

    def set_default_label(self):
        coord_system = self.get_grid_helper().get_display_coord_system()
        if coord_system is None:
            coord_system = self.get_grid_helper()._wcsgrid_orig_coord_system

        if coord_system == "fk5":
            xlabel=r"$\alpha_{2000}$"
            ylabel=r"$\delta_{2000}$"
        elif coord_system == "fk4":
            xlabel=r"$\alpha_{1950}$"
            ylabel=r"$\delta_{1950}$"
        elif coord_system == "gal":
            xlabel=r"$l$"
            ylabel=r"$b$"
        else:
            xlabel=r""
            ylabel=r""

        self.axis["left"].label.set_text(ylabel)
        self.axis["right"].label.set_text(ylabel)
        self.axis["bottom"].label.set_text(xlabel)
        self.axis["top"].label.set_text(xlabel)


    def swap_tick_coord(self):
        for axis in self.axis.values():
            gh = axis.get_helper()
            if isinstance(gh, FixedAxisArtistHelper):
                gh.change_tick_coord()

        label=self.axis["left"].label.get_text()
        self.axis["left"].label.set_text(self.axis["bottom"].label.get_text())
        self.axis["bottom"].label.set_text(label)

        label=self.axis["right"].label.get_text()
        self.axis["right"].label.set_text(self.axis["top"].label.get_text())
        self.axis["top"].label.set_text(label)

        self.get_grid_helper().invalidate()
        

SubplotWcs = maxes.subplot_class_factory(AxesWcs)



def test1():
    import pyfits, pywcs
    import matplotlib.pyplot as plt
    import axes_wcs
    fig = plt.figure(1)
    fig.clf()
    fname = "../doc/figures/data/lmc.fits"
    f = pyfits.open(fname)
    d, h = f[0].data, f[0].header

    global wcs
    wcs = pywcs.WCS(h)
    global grid_helper
    grid_helper = GridHelperWcs(wcs)
    global ax1
    ax1 = axes_wcs.SubplotWcs(fig, 1, 1, 1, grid_helper=grid_helper)

    fig.add_subplot(ax1)

    ax1.set_aspect(1.)
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 100)

    ax1.grid(True)

    plt.draw()


def test2():
    import pyfits, pywcs
    import matplotlib.pyplot as plt
    import axes_wcs

    fig = plt.figure(1)
    fig.clf()
    fname = "../doc/figures/data/lmc.fits"
    f = pyfits.open(fname)
    d, h = f[0].data, f[0].header

    ax1 = axes_wcs.SubplotWcs(fig, 1, 1, 1, header=h)

    fig.add_subplot(ax1)
    ax1.set_aspect(1.)

    ax1.grid(True)
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 100)

    #ax1.set_display_coord_system("gal")
    grid_helper_gal = ax1["gal"].get_grid_helper()
    grid_helper_gal.update_wcsgrid_params(coord_format=("dms","dms"))
    ax1.axis["b=-35"] = grid_helper_gal.new_floating_axis(1, -35, axes=ax1)
    ax1.axis["b=-35"].label.set_text("b = -35")

    plt.draw()


def test21():
    import pyfits, pywcs
    import matplotlib.pyplot as plt
    #import pywcsgrid2.axes_wcs as axes_wcs
    import axes_wcs
    #reload(axes_wcs)
    fig = plt.figure(1)
    fig.clf()
    #fname = "/Users/jjlee/local/src/astropy/pywcsgrid_all/test/i013b4h0.fit"
    fname = "../doc/figures/data/lmc.fits"
    f = pyfits.open(fname)
    d, h = f[0].data, f[0].header
    #d.shape = d.shape[-2:] # this particular image has a shape of [1, 500, 500]

    global wcs
    wcs = pywcs.WCS(h)
    #ax1 = axes_wcs.SubplotWcs(fig, 1, 1, 1, wcs=wcs)
    global grid_helper
    grid_helper = GridHelperWcs(wcs)
    ax1 = axes_wcs.SubplotWcs(fig, 1, 1, 1, grid_helper=grid_helper)

    fig.add_subplot(ax1)

    ax1.set_aspect(1.)
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 100)

    ax1.grid(True)

    #ax1.axis["top"].set_visible(False)
    #ax1.axis["right"].set_visible(False)

    global axx
    axx = ax1["fk5"].axis["top"]
    ax1["gal"].axis["top"].set_visible(True)
    #ax1["gal"].axis["top"].major_ticklabels.set_visible(True)
    ax1["gal"].axis["right"].set_visible(True)
    #ax1["gal"].axis["right"].major_ticklabels.set_visible(True)
    ax1["gal"].gridlines.set_visible(True)
    ax1["gal"].gridlines.set_color("r")

    #ax1["fk5"]._grid_helper = ax1.get_grid_helper()
    #new_fixed_axis = ax1["fk5"].get_grid_helper().new_fixed_axis
    #ax1["fk5"].axis["top2"] = new_fixed_axis(loc="top")
    #ax1["fk5"].axis["right"].set_visible(True)

    plt.draw()



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    #test1()
    test2()
    #test21()
    #test3()
    plt.show()
    #select_step(21.2, 33.3, 5)
    #select_step(20+21.2/60., 21+33.3/60., 5)
    #select_step(20.5+21.2/3600., 20.5+33.3/3600., 5)
    #select_step(20+21.2/60., 20+53.3/60., 5)
    #test_grider()

