
import matplotlib.axes as maxes
from matplotlib.cbook import is_string_like

import numpy as np
from pywcsgrid2.wcs_transforms import WcsSky2PixelTransform, WcsPixel2SkyTransform,\
     WcsSky2SkyTransform

from matplotlib.transforms import Affine2D

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
from wcs_helper import get_kapteyn_projection, coord_system_guess
#from pywcsgrid2.wcs_helper import coord_system_guess

from mpl_toolkits.axes_grid.parasite_axes import HostAxes, ParasiteAxesAuxTrans
import mpl_toolkits.axes_grid.grid_finder as grid_finder

GridFinderBase = grid_finder.GridFinderBase



from mpl_toolkits.axes_grid.angle_helper import ExtremeFinderCycle
from mpl_toolkits.axes_grid.angle_helper import LocatorDMS, LocatorHMS, FormatterDMS, FormatterHMS

from mpl_toolkits.axes_grid.grid_helper_curvelinear \
     import GridHelperCurveLinear, FixedAxisArtistHelper

from mpl_toolkits.axes_grid.anchored_artists import AnchoredEllipse, \
     AnchoredText, AnchoredSizeBar

from aux_artists import AnchoredCompass


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

class FormatterDMSDelta(object):
    def __call__(self, direction, factor, values):
        if len(values) == 0:
            return []
        if factor == 1:
            return ["$%d^{\circ}$" % (v,) for v in values]
        elif factor == 60:
            return ["$%d^{\prime}$" % (v,) for v in values]
        elif factor == 3600:
            return ["$%d^{\prime\prime}$" % (v,) for v in values]
        else: # factor > 3600.
            return ["$%s$" % (str(v),) for v in values]



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
            ctype1, ctype2 = self.projection.ctypes[:2]
            equinox = self.projection.equinox
            coord_guess = coord_system_guess(ctype1, ctype2, equinox)

            #coord_guess = coord_system_guess(wcs.wcs.ctype[0],
            #                                 wcs.wcs.ctype[1],
            #                                 wcs.wcs.equinox)
            if coord_guess is None:
                raise ValueError("Unknown coordinate system, %s, %s, equinox=%.1f" \
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

        #self._center_pixel = None
        self._center_world = None
        self._delta_trans = Affine2D()
        self._wcs_trans = None

        GridHelperCurveLinear.__init__(self, _wcs_trans,
                                       extreme_finder=ExtremeFinderCycle(20,20),
                                       grid_locator1=LocatorHMS(4),
                                       grid_locator2=LocatorDMS(4),
                                       tick_formatter1=FormatterHMS(),
                                       tick_formatter2=FormatterDMS())

#     def set_ticklabel_mode(self, m):
#         if m is None or m in ["center"]:
#             self._mode = m
#             #if self._center_pixel is None:
#             #    self.set_center_pixel(0, 0)
#             print "Center"
#         else:
#             raise ValueError("Unknown mode : " + str(m))

    def get_mode(self):
        return self._mode

#     def set_center_world(self, lon, lat):
#         x, y = self.projection.topixel(([lon], [lat]))
#         self._center_pixel = x[0]-1, y[0]-1
#         self._center_world = lon, lat

#     def tr(self, lon, lat):
#         if self.get_mode() == "center":
#             lon0, lat0 = self._center_world
#             lon, lat = lon + lon0, lat + lat0
#             x, y = self.projection.topixel((lon, lat))
#             x0, y0 = self._center_pixel
#             return (x-1)-x0, (y-1)-y0
#         else:
#             return x-1, y-1
#         #ll = np.concatenate((lon[:,np.newaxis], lat[:,np.newaxis]), 1)
#         #origin=0
#         #xy = self.wcs.wcs.s2p(ll, origin)["pixcrd"]
#         #return xy[:,0], xy[:,1]

#     def inv_tr(self, x, y):
#         if self.get_mode() == "center":
#             x0, y0 = self._center_pixel
#             x, y = x+x0, y+y0
#             lon, lat = self.projection.toworld((x+1, y+1))
#             lon0, lat0 = self._center_world
#             lon, lat = lon - lon0, lat - lat0
#             return lon, lat
#         else:
#             return self.projection.toworld((x+1, y+1))
#         #xy = np.concatenate((x[:,np.newaxis], y[:,np.newaxis]),
#         #                    1)
#         #origin=0
#         #ll = self.wcs.wcs.p2s(xy, origin)['world']
#         #return ll[:,0], ll[:,1]


    def _set_center_pixel(self, x, y):
        #self._center_pixel = x, y
        lon, lat = self.projection.toworld(([x+1], [y+1]))
        self._set_center_world(lon[0], lat[0])

    def _set_center_world(self, lon, lat):
        self._center_world = lon, lat
        cos_lat = np.cos(lat/180.*np.pi)
        self._delta_trans.clear().scale(-1./cos_lat,1.).translate(lon, lat)

    def set_ticklabel_mode(self, mode, **kwargs):
        if mode not in ["normal", "delta"]:
            raise ValueError("unknown mode : %s" % (mode))

        if mode == "delta":
            lon_lat = kwargs.pop("center_world", None)
            xy = kwargs.pop("center_pixel", None)
            if lon_lat and xy:
                raise ValueError("only one of center_world of center_pixel need to be set")
            if lon_lat:
                lon, lat = lon_lat
                self._set_center_world(lon, lat)
            elif xy:
                x, y = xy
                self._set_center_pixel(x, y)
            self.update_wcsgrid_params(coord_format=("delta", "delta"))
        else:
            self._delta_trans.clear()
            self.update_wcsgrid_params(coord_format=("hms", "dms"))

        if self._wcs_trans is not None:
            self.update_grid_finder(aux_trans=self._delta_trans+self._wcs_trans)
        else:
            _wcs_trans = self.get_wcs_trans(self.projection, None)
            self.update_grid_finder(aux_trans=self._delta_trans+_wcs_trans)


    def update_wcsgrid_params(self, **kwargs):
        """
        coord_format=("hms", "dms"),
        label_density=(4, 4),
        """

        self._wcsgrid_params.update(**kwargs)

        f1, f2 = self._wcsgrid_params["coord_format"]
        locator_selector = dict(hms=LocatorHMS,
                                dms=LocatorDMS,
                                delta=LocatorDMS).__getitem__
        formatter_selector = dict(hms=FormatterHMS,
                                  dms=FormatterDMS,
                                  delta=FormatterDMSDelta).__getitem__

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
        self._wcs_trans = self.get_wcs_trans(self.projection, coord_system)

        self.update_grid_finder(aux_trans=self._delta_trans+self._wcs_trans)


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
        """
        Creates an axes for displaying FITS images. The axes uses the
        FITS header information for ticks and grids appropriate for
        sky coordinate.

        The arguments for this axes is same as mpl's Axes, except that
        either *header* or *grid_helper* (but not both) keyword
        arguments must be provided.

        Necessary Keyword arguments:

          *header*: pyfits.Header instance
            header of the fits file that will be used for the axes
            a new pywcsgrid2.GridHelper instance is created using the header.

        or:

          *grid_helper*: pywcsgrid2.GridHelper instance


        """
        if "grid_helper" not in kw:
            header = kw.pop("header", None)
            wcs = kw.pop("wcs", None)
            if (header is not None) and (wcs is None):
                self.projection = get_kapteyn_projection(header)
            elif (header is None) and (wcs is not None):
                self.projection = get_kapteyn_projection(wcs)
            else:
                raise ValueError("wcs")

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
        """
        adjust the controlling option for ticks and grid.

        The option should be given as keyword arguments.


        Keyword arguments:

          *coord_format*: a tuple of two strings describing how ticks
            are located and labeld.  Currently,"hms", "dms" are only
            available options. The default is ("hms", "dms"), i.e.,
            (hour, minute, second) for the first coordinate and (degree,
            minute, second) for the second coordinate. For Galactic
            coordinate, ("dms", "dms") should be more sensible. The default
            is ("hms", "dms").


          *label_density*: approxmate number of ticks for each coordinates.
             Currently, this is the only parameter to control the
             tick locations. The default is (4, 4).

        """
        self.get_grid_helper().update_wcsgrid_params(**ka)


    def set_display_coord_system(self, c):
        """
        set the coordinate system to be displayed for the axes.

        Accept:
          "fk5"|"fk4"|"gal"
        """
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
        """
        Swap the coordinates that will be displayed for each axis.
        For example, by default, x axis shows ticklabels for the
        first coordinate(e.g., RA) and the y axis shows the second
        coordinate. swap_tick_coord makes the x-axis shows the
        second coordinate and y-axis the second. This would be only
        useful when you plot highly curved coordinate system.
        """
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


    def add_beam_size(self, major_pixel, minor_pixel, angle, loc,
                      frameon=False, borderpad=.8, patch_props=None,
                      **kwargs):
        """
        Add a ellipse patch to the axes with given sizes. This is
        intended to display beam size (or PSF size) of the observed
        image.  The ellipse is located according to the *loc*
        parameter, which should be the location code as in the legend.

        Keyword arguments:

          *major_pixel* : major axis size in pixel (in the image coordinate)

          *minor_pixel* : minor axis size in pixel

          *angle* : angle of rotation for the ellipse

          *patch_props* : a dictionary of properties to be passed for
             Patch creation, e.g. facecolor, alpha, etc.

        See mpl_toolkits.axes_grid.anchored_artists for orther options.
        """
        # Beam size
        # (major, minor) = 3, 4 in pixel, angle=20
        ae = AnchoredEllipse(self.transData, major_pixel, minor_pixel,
                             angle, loc=loc,
                             borderpad=borderpad, prop=None,
                             frameon=frameon, **kwargs)
        if patch_props is None:
            ae.ellipse.set(fc="none", ec="k") #, hatch="/")
        else:
            ae.ellipse.set(**patch_props)

        self.add_artist(ae)
        return ae


    def add_inner_title(self, title, loc, **kwargs):
        """
        Add a text at the inside corner of the axes.
        It simply uses mpl_toolkits.axes_grid.anchored_artists.AnchoredText.
        See mpl_toolkits.axes_grid.anchored_artists.AnchoredText
        for more details
        """
        # Figure title
        at = AnchoredText(title, loc=loc, **kwargs)
        self.add_artist(at)
        return at


    def add_compass(self, loc, coord="fk5", arrow_length=0.15,
                    txt1="E", txt2="N",
                    delta_a1=0, delta_a2=0, **kwargs):
        """
        Add two arros with appropriate labels showing the increasing
        direction of the coordinate. The default label is "E" and "N".


        Keyword arguments:

          *arrow_length* : length of the arrows in a fraction of axes size.

          *loc* : the location code as in the legend

          *coord* : the coordinate name. default is "fk5"

          *txt*, *txt2* : labels assocated with the arrows. Defaults are "E" & "N"

          delta_a1, delta_a2 : an optional angles (in degree) to be added to
            automatically determined directions.
        """
        # compass
        ac = AnchoredCompass(self, self[coord].transAux, loc=loc,
                             txt1=txt1, txt2=txt2,
                             delta_a1=delta_a1, delta_a2=delta_a2,
                             **kwargs)
        self.add_artist(ac)
        return ac


    def add_size_bar(self, length_pixel, label, loc, sep=5,
                     borderpad=0.8, frameon=False, **kwargs):
        """
        Add a horizontal line with a label underneath. Intended to
        display the angular size.

        Keyword arguments:

          *length_pixel* : length of the line in pixel in the image coordinate.

          *label* : label

          *loc* : the location code as in the legend

        """

        asb =  AnchoredSizeBar(self.transData,
                               length_pixel,
                               label,
                               loc=loc,
                               borderpad=borderpad, sep=sep,
                               frameon=frameon)
        self.add_artist(asb)

        return asb


SubplotWcs = maxes.subplot_class_factory(AxesWcs)




