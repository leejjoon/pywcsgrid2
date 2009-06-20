
import matplotlib.axes as maxes
from matplotlib.cbook import is_string_like
from matplotlib.path import Path

import numpy as np
from wcs_transforms import WcsSky2PixelTransform, WcsPixel2SkyTransform,\
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
from  mpl_toolkits.axes_grid.axislines import \
     AxisArtistHelper, GridHelperBase, AxisArtist, GridlinesCollection

import pywcs

from mpl_toolkits.axes_grid.parasite_axes import HostAxes, ParasiteAxesAuxTrans
import mpl_toolkits.axes_grid.grid_finder as grid_finder

GridFinderMplTransform = grid_finder.GridFinderMplTransform
GridFinderBase = grid_finder.GridFinderBase

from kapteyn_helper import coord_system_guess, sky2sky, is_equal_coord_sys
from kapteyn_helper import coord_system as _coord_sys_dict


from mpl_toolkits.axes_grid.angle_helper import ExtremeFinderCycle
from mpl_toolkits.axes_grid.angle_helper import LocatorDMS, LocatorHMS, FormatterDMS, FormatterHMS

from mpl_toolkits.axes_grid.grid_helper_curvelinear import FixedAxisArtistHelper

from mpl_toolkits.axes_grid.angle_helper import select_step24, select_step360


class GridFinderWcsFull(GridFinderBase):

    def __init__(self, wcs, tran,
                 extreme_finder,
                 grid_locator1,
                 grid_locator2,
                 tick_formatter1,
                 tick_formatter2):
        """
        transform : transfrom from the image coordinate (which will be
        the transData of the axes to the world coordinate.
        locator1, locator2 : grid locator for 1st and 2nd axis.
        """
        self.wcs = wcs
        self.lon_mode = ""
        self.lat_mode = ""


        if tran is not None:
            self.tran = tran
            self.tran_inv = tran.inverted()
        else:
            self.tran = self.tran_inv = None


        super(GridFinderWcsFull, self).__init__( \
                 extreme_finder,
                 grid_locator1,
                 grid_locator2,
                 tick_formatter1,
                 tick_formatter2)


    def get_grid_info(self,
                      x1, y1, x2, y2):

        return super(GridFinderWcsFull,self).get_grid_info( \
            self.transform_xy, self.inv_transform_xy,
            x1, y1, x2, y2)


    def transform_xy(self, x, y):
        lon, lat = self.wcs.wcs_pix2sky(x, y, 0)
        if self.tran is not None:
            lon, lat = self.tran(lon, lat)
        return lon, lat

    def inv_transform_xy(self, lon, lat):
        if self.tran_inv is not None:
            lon, lat = self.tran_inv(lon, lat)
        x, y = self.wcs.wcs_sky2pix(lon, lat, 0)
        return x, y



class GridFinderWcs(GridFinderWcsFull):

    def __init__(self, wcs, tran,
                 grid_mode=None):
        """
        transform : transfrom from the image coordinate (which will be
        the transData of the axes to the world coordinate.
        locator1, locator2 : grid locator for 1st and 2nd axis.
        """
        super(GridFinderWcs, self).__init__(wcs, tran,
                                            ExtremeFinderCycle(20,20),
                                            LocatorHMS(4),
                                            LocatorDMS(4),
                                            FormatterHMS(),
                                            FormatterDMS())






class GridHelperWcsBase(GridHelperBase):

    def _update(self, x1, x2, y1, y2):
        "bbox in 0-based image coordinates"
        # update wcsgrid

        if self._force_update is False and \
               self._old_values == (x1, x2, y1, y2,
                                    self.get_wcsgrid_params()):
            return

        if self._wcsgrid_display_coord_system is None:
            src_system = None
            tgt_system = None
        else:
            if self._wcsgrid_display_coord_system == self._wcsgrid_orig_coord_system:
                src_system = None
                tgt_system = None
            else:
                src_system = self._wcsgrid_orig_coord_system
                tgt_system = self._wcsgrid_display_coord_system


        self.update_wcsgrid_params(src_system=src_system,
                                   dest_system=tgt_system)


        self._update_grid(x1, y1, x2, y2,
                          **self._wcsgrid_params)

        self._old_values = (x1, x2, y1, y2, self.get_wcsgrid_params().copy())

        self._force_update = False


    def new_fixed_axis(self, loc,
                       nth_coord=None, passthrough_point=None,
                       tick_direction="in",
                       label_direction=None,
                       offset=None,
                       axes=None
                       ):

        _helper = FixedAxisArtistHelper(self, loc,
                                        nth_coord, nth_coord_ticks=None)

        axisline = AxisArtist(axes, _helper)

        return axisline



    # To be implemented in the derived class

    def get_wcsgrid_params(self):
        pass

    def _update_grid(self, x1, y1, x2, y2,
                     **wcsgrid_params):
        pass

    def get_gridlines(self):
        pass

    def get_tick_iterator(self, nth_coord, axis_side, minor=False):
        pass





class GridHelperWcs(GridHelperWcsBase):

    #def __init__(self, header):
    def __init__(self, wcs):

        super(GridHelperWcs, self).__init__()

        #self._header = header
        #wcs = pywcs.WCS(header)
        self._wcs = wcs
        self.wcsgrid = None
        self._old_values = None
        self._wcsgrid_params = dict(coord_format="GE",
                                    label_density=(6, 6),
                                    grid1=[], grid2=[])

        coord_guess = coord_system_guess(wcs.wcs.ctype[0],
                                         wcs.wcs.ctype[1],
                                         wcs.wcs.equinox)
        if coord_guess is None:
            raise ValueError("Unknown coordinate system, %, %s, equinox=%.1f" \
                             % (wcs.wcs.ctype[0], wcs.wcs.ctype[1],
                                wcs.wcs.equinox))
        else:
            self._wcsgrid_orig_coord_system = coord_guess

        self.set_display_coord_system(None)


        self.lon_locator = LocatorHMS(4)
        self.lat_locator = LocatorDMS(4)
        self.lon_formatter, self.lat_formatter = FormatterHMS, FormatterDMS


    def update_wcsgrid_params(self, **ka):
        """
        coord_format="GE",
        label_density=(6, 6),
        grid1=[], grid2=[]):
        """

        for k in ["coord_format", "label_density", "grid1", "grid2", "src_system", "dest_system"]:
            if k in ka:
                self._wcsgrid_params[k] = ka.pop(k)

        if len(ka) > 0:
            raise ValueError("keyword name should be one of coord_format, label_densit, grid1, grid2")


    def get_wcsgrid_params(self):
        return self._wcsgrid_params

    def _update_grid(self, x1, y1, x2, y2,
                     **wcsgrid_params):

        if is_equal_coord_sys(self._wcsgrid_orig_coord_system,
                              self._wcsgrid_display_coord_system):
            tran = None
        else:
            tran = sky2sky(self._wcsgrid_orig_coord_system,
                           self._wcsgrid_display_coord_system)

        grider = GridFinderWcs(self._wcs, tran)
        self.wcsgrid = grider.get_grid_info(x1, y1, x2, y2)


    def get_gridlines(self):
        grid_lines = []
        for gl in self.wcsgrid["lat"]["lines"]:
            grid_lines.extend(gl)
        for gl in self.wcsgrid["lon"]["lines"]:
            grid_lines.extend(gl)

        return grid_lines


    def get_tick_iterator(self, nth_coord, axis_side, minor=False):

        axisnr = dict(left=0, bottom=1, right=2, top=3)[axis_side]
        angle = [0, 90, 180, 270][axisnr]
        #side_select = ["left", "bottom", "right", "top"].__getitem__
        lon_or_lat = ["lon", "lat"][nth_coord]
        if not minor: # major ticks
            def f():
                for (xy, a), l in zip(self.wcsgrid[lon_or_lat]["tick_locs"][axis_side],
                                    self.wcsgrid[lon_or_lat]["tick_labels"][axis_side]):
                    yield xy, a, l
        else:
            def f():
                for (xy, a), l in zip(self.wcsgrid[lon_or_lat]["tick_locs"][axis_side],
                                    self.wcsgrid[lon_or_lat]["tick_labels"][axis_side]):
                    yield xy, a, ""
                #for xy, a, l in self.wcsgrid[lon_or_lat]["ticks"][axis_side]:
                #    yield xy, a, ""

        return f()


    def set_display_coord_system(self, coord_system):
        if self._wcsgrid_orig_coord_system is None:
            raise ValueError("_pywcsgrid_orig_coord_system is not set")

        if coord_system is None:
            coord_system = self._wcsgrid_orig_coord_system

        self._wcsgrid_display_coord_system = coord_system
        self.invalidate()

    def set_display_coord_system2(self, coord_system):
        if self._wcsgrid_orig_coord_system is None:
            raise ValueError("_pywcsgrid_orig_coord_system is not set")

        if coord_system is None:
            coord_system = self._wcsgrid_orig_coord_system

        if is_string_like(coord_system):
            coord_system = _coord_sys_dict[coord_system.lower()]

        self._wcsgrid_display_coord_system = coord_system
        self.invalidate()

    def get_display_coord_system(self):
        return self._wcsgrid_display_coord_system



class ParasiteAxesSky(ParasiteAxesAuxTrans):
    def __init__(self, parent, src_coord, *kl, **kwargs):

        dest_coord = parent.get_grid_helper()._wcsgrid_orig_coord_system
        tr = WcsSky2SkyTransform(src_coord, dest_coord) + \
             parent._wcsgrid_wcsaxes[0].transAux


        grid_helper = GridHelperWcs(parent._wcs)
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
        #return super(ParasiteAxesWcs, self).imshow(Om, extent=oe, **kwargs)

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
            wcs = kw.pop("wcs", None)
            if wcs is None:
                raise ValueError("wcs")
            else:
                self._wcs = wcs

            grid_helper = GridHelperWcs(self._wcs)
            kw["grid_helper"] = grid_helper
        else:
            self._wcs = kw["grid_helper"]._wcs

        super(AxesWcs, self).__init__(*kl, **kw)

        self._init_parasites()

        self.set_default_label()


    def _init_parasites(self):
        ax = ParasiteAxesAuxTrans(self,
                                  WcsSky2PixelTransform(self._wcs),
                                  viewlim_mode="equal")
        self._wcsgrid_wcsaxes = {0:ax}
        self.parasites.append(ax)


    def register_wcs(self, key, wcs):
        pass

    def __getitem__(self, key):

        # check if key is a valid coord_sys instance
        if key == 0: # or isinstance(key, _coord_sys_dict) :
            pass
        elif is_string_like(key):
            pass
            #key = _coord_sys_dict[key.lower()]
        elif isinstance(key, pywcs.WCS):
            pass
        else:
            raise ValueError("invalide key : %s" % repr(key))

        if key not in self._wcsgrid_wcsaxes:
            if self.get_grid_helper()._wcsgrid_orig_coord_system == key:
                self._wcsgrid_wcsaxes[key] = self._wcsgrid_wcsaxes[0]
            else:
                orig_coord = self.get_grid_helper()._wcsgrid_orig_coord_system

                if isinstance(key, pywcs.WCS):
                    ax = ParasiteAxesWcs(self, key)
                else:
                    ax = ParasiteAxesSky(self, key)

                self._wcsgrid_wcsaxes[key] = ax
                self.parasites.append(ax)


        return self._wcsgrid_wcsaxes[key]

    def set_display_coord_system(self, c):
        self.get_grid_helper().set_display_coord_system(c)

        self.axis["bottom"].get_helper().change_tick_coord(0)
        self.axis["top"].get_helper().change_tick_coord(0)
        self.axis["left"].get_helper().change_tick_coord(1)
        self.axis["right"].get_helper().change_tick_coord(1)

        self.set_default_label()

    def set_default_label(self):
        coord_system = self.get_grid_helper().get_display_coord_system()
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
            axis.get_helper().change_tick_coord()

        label=self.axis["left"].label.get_text()
        self.axis["left"].label.set_text(self.axis["bottom"].label.get_text())
        self.axis["bottom"].label.set_text(label)

        label=self.axis["right"].label.get_text()
        self.axis["right"].label.set_text(self.axis["top"].label.get_text())
        self.axis["top"].label.set_text(label)



SubplotWcs = maxes.subplot_class_factory(AxesWcs)




def test21():
    import pyfits, pywcs
    import matplotlib.pyplot as plt
    import axes_wcs
    reload(axes_wcs)
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
    #test2()
    test21()
    #test3()
    plt.show()
    #select_step(21.2, 33.3, 5)
    #select_step(20+21.2/60., 21+33.3/60., 5)
    #select_step(20.5+21.2/3600., 20.5+33.3/3600., 5)
    #select_step(20+21.2/60., 20+53.3/60., 5)
    #test_grider()

