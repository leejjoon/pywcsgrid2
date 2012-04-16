import warnings

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


#import pywcs
from wcs_helper import get_kapteyn_projection, coord_system_guess
#from pywcsgrid2.wcs_helper import coord_system_guess

from  mpl_toolkits.axisartist import GridlinesCollection
import mpl_toolkits.axisartist.grid_finder as grid_finder

GridFinderBase = grid_finder.GridFinderBase

from mpl_toolkits.axisartist.angle_helper import ExtremeFinderCycle

from locator_formatter import LocatorDMS, LocatorHMS, FixedLocator, \
     FixedFormatter, MaxNLocator, FormatterPrettyPrint



from mpl_toolkits.axisartist.angle_helper import FormatterDMS, FormatterHMS

from mpl_toolkits.axisartist.grid_helper_curvelinear \
     import FixedAxisArtistHelper

import mpl_toolkits.axisartist.grid_helper_curvelinear as grid_helper_curvelinear


# implment locator_params that is not ioncluded in matplotlib 1.0
class GridHelperCurveLinear(grid_helper_curvelinear.GridHelperCurveLinear):
    def locator_params(self, axis='both', **kwargs):
        """
        Convenience method for controlling tick locators.

        Keyword arguments:

        *axis*
            ['x' | 'y' | 'both']  Axis on which to operate;
            default is 'both'.

        Remaining keyword arguments are passed to directly to
        set_params method of each locators.

        Typically one might want to reduce the maximum number
        of ticks, for example::

            ax.locator_params(nbins=4)

        """

        if kwargs.has_key("tight"):
            warnings.warn("axislines ignores *tight* parameter")

        _x = axis in ['x', 'both']
        _y = axis in ['y', 'both']

        if _x:
            self.grid_finder.grid_locator1.set_params(**kwargs)
        if _y:
            self.grid_finder.grid_locator2.set_params(**kwargs)

        self.invalidate()



from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse, \
     AnchoredText, AnchoredSizeBar

from aux_artists import AnchoredCompass

import mpl_toolkits.axisartist as axisartist
from mpl_toolkits.axisartist import ParasiteAxesAuxTrans
import mpl_toolkits.axisartist.floating_axes as floating_axes

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

import matplotlib.ticker as mticker

class FormatterDMSDelta(object):
    def __init__(self, usetex=True):
        self._fmt = mticker.ScalarFormatter()
        self._fmt._usetex = usetex
        self._fmt.create_dummy_axis()
        self._ignore_factor = True

    def _get_pretty_fromat(self, direction, factor, values):
        if not self._ignore_factor:
            if factor is None:
                factor = 1.
            values = [v/factor for v in values]
        #values = [v for v in values]
        self._fmt.set_locs(values)
        return [self._fmt(v) for v in values]


    def __call__(self, direction, factor, values):
        values = self._get_pretty_fromat(direction, factor, values)
        if len(values) == 0:
            return []
        if factor == 1 or factor is None:
            return ["%s$^{\circ}$" % (v,) for v in values]
        elif factor == 60:
            return ["%s$^{\prime}$" % (v,) for v in values]
        elif factor == 3600:
            return ["%s$^{\prime\prime}$" % (v,) for v in values]
        elif factor > 3600.:
            return ["%s" % (v,) for v in values]
        else: # factor > 3600.
            return ["%s" % (v,) for v in values]


#class SimpleScaleTranslateAffine2D(Affine2D):
class SimpleScaleTranslateAffine2D(Affine2D):
    """
    A mutable 2D affine transformation.
    """

    def __init__(self):
        """
        Initialize an Affine transform from a 3x3 numpy float array::

          a c e
          b d f
          0 0 1

        If *matrix* is None, initialize with the identity transform.
        """
        Affine2D.__init__(self)
        matrix = np.identity(3)
        self._mtx = matrix
        self._invalid = 0

    def __repr__(self):
        return "Affine2D(%s)" % repr(self._mtx)
    __str__ = __repr__

    def __cmp__(self, other):
        if (isinstance(other, Affine2D) and
            (self.get_matrix() == other.get_matrix()).all()):
            return 0
        return -1

    @staticmethod
    def from_values(a, d, e, f):
        """
        (staticmethod) Create a new Affine2D instance from the given
        values::

          a c e
          b d f
          0 0 1
        """
        b, c = 0, 0
        return Affine2D(
            np.array([a, c, e, b, d, f, 0.0, 0.0, 1.0], np.float_)
            .reshape((3,3)))

    def get_matrix(self):
        """
        Get the underlying transformation matrix as a 3x3 numpy array::

          a c e
          b d f
          0 0 1
        """
        self._invalid = 0
        return self._mtx

    def set_matrix(self, mtx):
        """
        Set the underlying transformation matrix from a 3x3 numpy array::

          a c e
          b d f
          0 0 1
        """
        self._mtx = mtx
        self.invalidate()

    def set(self, other):
        """
        Set this transformation from the frozen copy of another
        :class:`Affine2DBase` object.
        """
        assert isinstance(other, Affine2D)
        self._mtx = other.get_matrix()
        self.invalidate()

    def clear(self):
        """
        Reset the underlying matrix to the identity transform.
        """
        self._mtx = np.identity(3)
        self.invalidate()
        return self

    def set_tx(self, tx):
        """
        set a translation x in place.

        """
        self._mtx[0,2] = tx
        self.invalidate()

    def set_ty(self, ty):
        """
        set a translation in place.

        """
        self._mtx[1,2] = ty
        self.invalidate()

    def set_sx(self, sx):
        """
        Adds a scale x in place.

        """

        self._mtx[0, 0] = sx
        self.invalidate()

    def set_sy(self, sy):
        """
        Adds a scale x in place.

        """

        self._mtx[1, 1] = sy
        self.invalidate()

    def reset_x(self):
        """
        reset x in place.

        """

        self.set_tx(0)
        self.set_sx(1)
        self.invalidate()

    def reset_y(self):
        """
        reset x in place.

        """

        self.set_ty(0)
        self.set_sy(1)
        self.invalidate()


    def _get_is_separable(self):
        mtx = self.get_matrix()
        return mtx[0, 1] == 0.0 and mtx[1, 0] == 0.0
    is_separable = property(_get_is_separable)



import re
_pattern_longitude = re.compile(r"^(RA|GLON)")
_pattern_latitude = re.compile(r"^(DEC|GLAT)")



class GridHelperWcsBase(object):

    # derived class class must define "_GRIDHELPER_CLASS"

    WCS_TRANS_COLLECTION = wcs_trans_collection

    @classmethod
    def get_wcs_trans(cls, wcs, coord):
        return cls.WCS_TRANS_COLLECTION.get_wcs_trans(wcs, coord)

    def _init_wcsgrid_params(self):
        self._wcsgrid_params = dict(coord_format=("hms", "dms"),
                                    label_density=(4, 4))

    def _init_projection(self, wcs, axis_nums=None):
        #self.header = header
        #wcs = pywcs.WCS(header)
        projection = get_kapteyn_projection(wcs)

        if axis_nums is None:
            axis_nums = [0,1]
        self.projection = projection.substitute(axis_nums, [0]*projection.naxis)

        #self.projection = projection

        self.wcsgrid = None
        self._old_values = None

        self._init_wcsgrid_params()
        self._wcs_trans = self.get_wcs_trans(self.projection, None)
        self._wcsgrid_orig_coord_system = None

    def __init__(self, wcs, axis_nums=None):
        self._init_projection(wcs, axis_nums)

        #self._center_world = None
        #self._delta_adef = 1, 1, 0, 0
        #self._delta_trans = SimpleScaleTranslateAffine2D()
        self._delta_trans = SimpleScaleTranslateAffine2D()

    def update_delta_trans(self, tx=None, ty=None, sx=None, sy=None):
        if tx is not None:
            self._delta_trans.set_tx(tx)
        if ty is not None:
            self._delta_trans.set_ty(ty)
        if sx is not None:
            self._delta_trans.set_sx(sx)
        if sy is not None:
            self._delta_trans.set_sy(sy)

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


    def _set_center_pixel_deprecated(self, x, y):
        #self._center_pixel = x, y
        lon, lat = self.projection.toworld(([x+1], [y+1]))
        cos_lat = np.cos(lat[0]/180.*np.pi)
        self._set_center_world(lon[0], lat[0],
                               lon_scale=1./cos_lat,
                               lat_scale=1.)

    def _set_center_world_deprecated(self, lon, lat, lon_scale, lat_scale):
        self._center_world_and_scale = lon, lat, lon_scale, lat_scale
        #cos_lat = np.cos(lat/180.*np.pi)
        self._delta_trans.clear().scale(lon_scale, lat_scale).translate(lon, lat)


    def set_ticklabel_mode_deprecated(self, mode, **kwargs):
        if mode not in ["normal", "delta"]:
            raise ValueError("unknown mode : %s" % (mode))

        if mode == "delta":
            lon_lat = kwargs.pop("center_world", None)
            xy = kwargs.pop("center_pixel", None)
            if lon_lat and xy:
                raise ValueError("only one of center_world of center_pixel need to be set")
            if lon_lat:
                lon, lat = lon_lat
                cos_lat = np.cos(lat/180.*np.pi)
                self._set_center_world(lon, lat,
                                       lon_scale=1./cos_lat,
                                       lat_scale=1.)
            elif xy:
                x, y = xy
                self._set_center_pixel_deprecated(x, y)
            else:
                raise RuntimeError("center_world of center_pixel argument is required")

            self.update_wcsgrid_params(coord_format=("delta", "delta"))
        else:
            self._delta_trans.clear()
            self._center_world = None
            self.update_wcsgrid_params(coord_format=("hms", "dms"))

        if self._wcs_trans is not None:
            self.update_grid_finder(aux_trans=self._delta_trans+self._wcs_trans)
        else:
            _wcs_trans = self.get_wcs_trans(self.projection, None)
            self.update_grid_finder(aux_trans=self._delta_trans+_wcs_trans)

    def update_wcsgrid_params(self, **kwargs):
        """
        This method is deprecated. Use set_ticklabel_type or locator_params.

        coord_format=("hms", "dms"),
        label_density=(4, 4),
        """

        warnings.warn("update_wcsgrid_params is deprecated. Use set_ticklabel_type or locator_params")

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



    def _get_default_locator_formatter(self, ctype_name):
        ctype_name = ctype_name.upper()
        if ctype_name.startswith("RA"):
            return LocatorHMS(4), FormatterHMS()
        elif ctype_name.startswith("DEC"):
            return LocatorDMS(4), FormatterDMS()
        elif ctype_name.startswith("GLON") or ctype_name.startswith("GLAT"):
            return LocatorDMS(4), FormatterDMS()
        else:
            return None, None


    def _get_default_extreme_finder_par(self, ctype_name):
        ctype_name = ctype_name.upper()
        if ctype_name.startswith("RA") or ctype_name.startswith("GLON"):
            cycle=360.0
            minmax=None
        elif ctype_name.startswith("DEC") or ctype_name.startswith("GLAT"):
            cycle=None
            minmax=None
        else:
            cycle=None
            minmax=None

        return cycle, minmax

    def _get_default_grid_helper_params(self, ctype1, ctype2,
                                        extreme_finder=None,
                                        grid_locator1=None,
                                        grid_locator2=None,
                                        tick_formatter1=None,
                                        tick_formatter2=None,
                                        ):
        if extreme_finder is None:
            lon_cycle, lon_minmax = self._get_default_extreme_finder_par(ctype1)
            lat_cycle, lat_minmax = self._get_default_extreme_finder_par(ctype2)

            extreme_finder= ExtremeFinderCycle(20, 20,
                                               lon_cycle=lon_cycle, lat_cycle=lat_cycle,
                                               lon_minmax=lon_minmax, lat_minmax=lat_minmax)

        _grid_locator1, _tick_formatter1 = self._get_default_locator_formatter(ctype1)
        if grid_locator1 is None:
            grid_locator1 = _grid_locator1
        if tick_formatter1 is None:
            tick_formatter1 = _tick_formatter1

        _grid_locator2, _tick_formatter2 = self._get_default_locator_formatter(ctype2)
        if grid_locator2 is None:
            grid_locator2 = _grid_locator2
        if tick_formatter2 is None:
            tick_formatter2 = _tick_formatter2

        return dict(extreme_finder=extreme_finder,
                    grid_locator1=grid_locator1,
                    grid_locator2=grid_locator2,
                    tick_formatter1=tick_formatter1,
                    tick_formatter2=tick_formatter2)


    def _set_ticklabel_type(self, labtyp, **labtyp_kwargs):

        nbins = labtyp_kwargs.pop("nbins", 4)
        locs = labtyp_kwargs.pop("locs", None)

        if labtyp == "manual":
            offset = labtyp_kwargs.pop("offset", 0)
            scale = labtyp_kwargs.pop("scale", 1)

            labels = labtyp_kwargs.pop("labels", None)

            if locs is None or labels is None:
                raise ValueError("manual labtyp requires locs and labels parameters")
            if len(locs) != len(labels):
                raise ValueError("locs and labels must be losts of same lengths.")


            locator = FixedLocator(locs)
            locator.set_factor(scale)

            formatter = FixedFormatter(labels)

            return offset, scale, locator, formatter

        elif labtyp in ["delta", "arcdeg", "arcmin", "arcsec"]:
            offset = labtyp_kwargs.pop("offset")
            scale = labtyp_kwargs.pop("scale", 1)

            if labtyp in ["delta"]:
                locator = LocatorDMS(nbins)
            elif labtyp in ["arcdeg", "arcmin", "arcsec"]:
                if locs is not None:
                    locator = FixedLocator(locs)
                else:
                    locator = MaxNLocator(nbins)

            if labtyp == "arcmin":
                locator.set_factor(60.)
            elif labtyp == "arcsec":
                locator.set_factor(3600.)

            if labtyp in ["delta"]:
                formatter = FormatterDMSDelta()
            else:
                formatter = FormatterPrettyPrint()

        else:
            if labtyp in ["absval", "absdeg"]:
                if locs is not None:
                    locator = FixedLocator(locs)
                else:
                    locator = MaxNLocator(nbins)

                if labtyp == "absval":
                    scale = labtyp_kwargs.pop("scale", 1)
                    locator.set_factor(scale)
                elif labtyp == "absdeg":
                    pass

                formatter = FormatterPrettyPrint()

            elif labtyp == "hms":
                locator = LocatorHMS(nbins)
                formatter = FormatterHMS()
            elif labtyp == "dms":
                locator = LocatorDMS(nbins)
                formatter = FormatterDMS()
#             else:
#                 #self.update_delta_trans(ty=0, sy=1)
#                 locator = MaxNLocator(nbins)
#                 formatter = FormatterPrettyPrint()
            else:
                raise ValueError("Unknown labeltype : %s" % (labtyp))

            offset, scale = 0, 1

        return offset, scale, locator, formatter


    def _check_scale_for_longitude(self, naxis, labtyp, labtyp_kwargs):
        if labtyp in ["delta", "arcdeg", "arcmin", "arcsec"]:
            ctype = self.projection.ctypes[naxis]
            if _pattern_longitude.match(ctype) and "scale" not in labtyp_kwargs:
                lat = labtyp_kwargs.pop("latitude", None)
                if lat is None:
                    raise RuntimeError("For ctype of longitude, either scale or the latitude parameter needs to be set")
                else:
                    scale = 1./np.cos(lat/180.*np.pi)
                    labtyp_kwargs["scale"] = scale


    def set_ticklabel1_type(self, labtyp, **labtyp_kwargs):
        """
        labtyp : there are two kinds of labtyp, offset-types and
        non-offset-type.  For offset-type labtyp, the label values are
        shifted from the original by the amount of the given
        *offset*. For either types, one can optionally specify the
        *scale* parameter and the scaled values (offset will be
        applied first if given) will be used. i.e.,

            world1 = (world0 - offset) * scale

        Followings are currently available labtyp.

         * non-offset : absval, absdeg, hms, dms
         * offset : delta, arcdeg, arcmin, arcsec


        By default, the tick locations will be chosen
        automatically. And *nbins* parameter controls the approximate
        number of ticks. When you want your ticks at the specified
        location, use the *locs* parameter, which takes a list of tick
        values. However, manual ticks are not possible for labtyp of
        hms, dms, delta. The *locs* parameter will be simply ignored.

        For the wcs coordinate, the default values are in the unit of
        degree.  And for the ctype of longitude, you may specify
        latitude (in degree) instead of the scale.

        .. plot:: figures/demo_labtyp.py

        """
        self._check_scale_for_longitude(0, labtyp, labtyp_kwargs)

        offset, scale, locator, formatter = \
                self._set_ticklabel_type(labtyp, **labtyp_kwargs)
        self.update_delta_trans(sx=scale, tx=offset)

        self.update_grid_finder(aux_trans=self._delta_trans+self._wcs_trans,
                                grid_locator1=locator,
                                tick_formatter1=formatter
                                )

    def set_ticklabel2_type(self, labtyp, **labtyp_kwargs):
        """
        see :meth:`~pywcsgrid2.axes_wcs.GridHelperWcsBase.set_ticklabel1_type`.
        """
        self._check_scale_for_longitude(1, labtyp, labtyp_kwargs)

        offset, scale, locator, formatter = \
                self._set_ticklabel_type(labtyp, **labtyp_kwargs)
        self.update_delta_trans(sy=scale, ty=offset)

        self.update_grid_finder(aux_trans=self._delta_trans+self._wcs_trans,
                                grid_locator2=locator,
                                tick_formatter2=formatter
                                )



    def _get_center_world_and_default_scale_factor(self, center_pixel):
        x, y = center_pixel
        center_world = self.projection.toworld(([x+1], [y+1]))
        clon, clat = center_world[0][0], center_world[1][0]

        ctype1, ctype2 = self.projection.ctypes[:2]

        if _pattern_latitude.match(ctype2):
            scale2 = 1
            if _pattern_longitude.match(ctype1):
                scale1 = 1./np.cos(clat/180.*np.pi)
            else:
                scale1 = None
        if _pattern_latitude.match(ctype1):
            scale1 = 1
            if _pattern_longitude.match(ctype2):
                scale2 = 1./np.cos(clon/180.*np.pi)
            else:
                scale2 = None

        return clon, clat, scale1, scale2



    def set_ticklabel_type(self, labtyp1, labtyp2,
                           center_pixel=None,
                           labtyp1_kwargs=None, labtyp2_kwargs=None,
                           **kwargs):
        """
        Set the label type of the x and y axis. Available options are

        =========  ================================================
         "hms"     label is in H M S.S (e.g. for RA)
         "dms"     label is in D M S.S (e.g. for DEC)
         "absval"  label is in actual wcs value
         "absdeg"  label is in degrees
         "delta"   label is in offsets (arcsec, arcmin, etc)
         "arcdeg"  label is in arcdegree offsets
         "arcmin"  label is in arcminute offsets
         "arcsec"  label is in arcsecond offsets
        =========  ================================================

        If labtyp2 is None, labtyp1 is used.

        For labtyp in ["delta", "arcdeg", "arcmin", "arcsec"],
        *center_pixel* parameter needs to be specified, from which
        the offsets will be calculated.

         center_pixel : in pixel coordinate (data coordinate in matplotlib)

        If None, the center position of the current view-limit will be used.

        *labtyp1_kwargs* and *labtyp2_kwargs* should be a dictionary of
        keyword arguments that are associated with *labtyp1* and
        *labtyp2*. Any other keyword arguments will be associated both with
        *labtyp1* and *labtyp2*. The available keywords arguments will
        depends on the value of *labtyp1* and *labtyp2*. Most common
        options are *nbins* and *locs*. *nbins* specify approximate number
        of automatically generated ticks. *locs* manually
        specify the tick locations.


        labtyp of individual axis can be set using
        :meth:`~pywcsgrid2.axes_wcs.GridHelperWcsBase.set_ticklabel1_type`
        and
        :meth:`~pywcsgrid2.axes_wcs.GridHelperWcsBase.set_ticklabel2_type`
        methods. See documentaion of these methods for more details.
        """

        if center_pixel is  None:
            tx, ty, sx, sy = None, None, None, None
        else:
            tx, ty, sx, sy = self._get_center_world_and_default_scale_factor(center_pixel)

        if labtyp1_kwargs is None:
            labtyp1_kwargs = dict(offset=tx, scale=sx)
        else:
            if "scale" not in labtyp1_kwargs:
                labtyp1_kwargs["scale"] = sx
            if "offset" not in labtyp1_kwargs:
                labtyp1_kwargs["offset"] = tx

        if labtyp2_kwargs is None:
            labtyp2_kwargs = dict(offset=ty, scale=sy)
        else:
            if "scale" not in labtyp2_kwargs:
                labtyp2_kwargs["scale"] = sy
            if "offset" not in labtyp2_kwargs:
                labtyp2_kwargs["offset"] = ty

        self.set_ticklabel1_type(labtyp1, **labtyp1_kwargs)
        self.set_ticklabel2_type(labtyp2, **labtyp2_kwargs)



class GridHelperWcsSkyBase(GridHelperWcsBase):
    _GRIDHELPER_CLASS=GridHelperCurveLinear

    def _init_projection(self, wcs, orig_coord, axis_nums=None):

        GridHelperWcsBase._init_projection(self, wcs, axis_nums=None)

        if orig_coord is None:
            ctype1, ctype2 = self.projection.ctypes[:2]
            equinox = self.projection.equinox
            coord_guess = coord_system_guess(ctype1, ctype2, equinox)

            #coord_guess = coord_system_guess(wcs.wcs.ctype[0],
            #                                 wcs.wcs.ctype[1],
            #                                 wcs.wcs.equinox)

            self._wcsgrid_orig_coord_system = coord_guess
            if coord_guess is None:

                print "Unknown coordinate system, %s, %s, equinox=%.1f" \
                                 % (ctype1, ctype2, equinox)
            #else:
            #    self._wcsgrid_orig_coord_system = coord_guess
        else:
            self._wcsgrid_orig_coord_system = orig_coord


        self._wcsgrid_display_coord_system = None


    def __init__(self, wcs, orig_coord=None, axis_nums=None):
        self._init_projection(wcs, orig_coord, axis_nums)

        self._center_world = None
        #self._delta_adef = 1, 1, 0, 0
        #self._delta_trans = SimpleScaleTranslateAffine2D()
        self._delta_trans = SimpleScaleTranslateAffine2D()


    def set_display_coord_system(self, coord_system):
        if self._wcsgrid_orig_coord_system is None:
            raise ValueError("_pywcsgrid_orig_coord_system is not set")

        self._wcsgrid_display_coord_system = coord_system
        self._wcs_trans = self.get_wcs_trans(self.projection, coord_system)

        if coord_system in ["fk5", "fk4"]:
            ctype1, ctype2 = "RA", "DEC"
        elif coord_system == "gal":
            ctype1, ctype2 = "GLON", "GLAT"
        else:
            ctype1, ctype2 = None, None

        grid_helper_params = self._get_default_grid_helper_params(ctype1, ctype2)

        self.update_grid_finder(aux_trans=self._delta_trans+self._wcs_trans,
                                grid_locator1=grid_helper_params["grid_locator1"],
                                grid_locator2=grid_helper_params["grid_locator2"],
                                tick_formatter1=grid_helper_params["tick_formatter1"],
                                tick_formatter2=grid_helper_params["tick_formatter2"]
                                )


    def get_display_coord_system(self):
        return self._wcsgrid_display_coord_system


class GridHelperWcsSimple(GridHelperWcsBase, GridHelperCurveLinear):
    _GRIDHELPER_CLASS=GridHelperCurveLinear

    def __init__(self, wcs,
                 extreme_finder=None,
                 grid_locator1=None,
                 grid_locator2=None,
                 tick_formatter1=None,
                 tick_formatter2=None,
                 axis_nums=None,
                 ):

        GridHelperWcsBase.__init__(self, wcs, axis_nums)

        ctype1, ctype2 = self.projection.ctypes[:2]

        grid_helper_params = self._get_default_grid_helper_params(ctype1, ctype2,
                                                                  extreme_finder=extreme_finder,
                                                                  grid_locator1=grid_locator1,
                                                                  grid_locator2=grid_locator2,
                                                                  tick_formatter1=tick_formatter1,
                                                                  tick_formatter2=tick_formatter2)

        GridHelperCurveLinear.__init__(self, self._wcs_trans,
                                       extreme_finder=grid_helper_params["extreme_finder"],
                                       grid_locator1=grid_helper_params["grid_locator1"],
                                       grid_locator2=grid_helper_params["grid_locator2"],
                                       tick_formatter1=grid_helper_params["tick_formatter1"],
                                       tick_formatter2=grid_helper_params["tick_formatter2"])



class GridHelperWcsSky(GridHelperWcsSkyBase, GridHelperCurveLinear):
    _GRIDHELPER_CLASS=GridHelperCurveLinear

    def __init__(self, wcs, orig_coord=None,
                 extreme_finder=None,
                 grid_locator1=None,
                 grid_locator2=None,
                 tick_formatter1=None,
                 tick_formatter2=None,
                 axis_nums=None,
                 ):

        GridHelperWcsSkyBase.__init__(self, wcs, orig_coord, axis_nums)

        # if extreme_finder is None:
        #     extreme_finder=ExtremeFinderCycle(20,20)
        # if grid_locator1 is None:
        #     grid_locator1=LocatorHMS(4)
        # if grid_locator2 is None:
        #     grid_locator2=LocatorDMS(4)
        # if tick_formatter1 is None:
        #     tick_formatter1=FormatterHMS()
        # if tick_formatter2 is None:
        #     tick_formatter2=FormatterDMS()

        #ctype1 = proj.ctypes[axis1]
        #ctype2 = proj.ctypes[axis2]

        if orig_coord in ["fk5", "fk4"]:
            ctype1, ctype2 = "RA", "DEC"
        elif orig_coord == "gal":
            ctype1, ctype2 = "GLON", "GLAT"
        elif orig_coord is None:
            ctype1, ctype2 = self.projection.ctypes[:2]
        else:
            raise ValueError("unknown coord type")


        grid_helper_params = self._get_default_grid_helper_params(ctype1, ctype2,
                                                                  extreme_finder=extreme_finder,
                                                                  grid_locator1=grid_locator1,
                                                                  grid_locator2=grid_locator2,
                                                                  tick_formatter1=tick_formatter1,
                                                                  tick_formatter2=tick_formatter2)

        GridHelperCurveLinear.__init__(self, self._wcs_trans,
                                       extreme_finder=grid_helper_params["extreme_finder"],
                                       grid_locator1=grid_helper_params["grid_locator1"],
                                       grid_locator2=grid_helper_params["grid_locator2"],
                                       tick_formatter1=grid_helper_params["tick_formatter1"],
                                       tick_formatter2=grid_helper_params["tick_formatter2"])

    #def _init_wcsgrid_params(self):
    #    self._wcsgrid_params = dict(coord_format=("hms", "dms"),
    #                                label_density=(4, 4))


GridHelperWcs = GridHelperWcsSky


class GridHelperWcsFloating(GridHelperWcsBase, floating_axes.GridHelperCurveLinear):
    def __init__(self, wcs, extremes,
                 orig_coord=None,
                 grid_locator1=None,
                 grid_locator2=None,
                 tick_formatter1=None,
                 tick_formatter2=None,
                 ):

        GridHelperWcsBase.__init__(self, wcs, orig_coord)


        if grid_locator1 is None:
            grid_locator1=LocatorHMS(4)
        if grid_locator2 is None:
            grid_locator2=LocatorDMS(4)
        if tick_formatter1 is None:
            tick_formatter1=FormatterHMS()
        if tick_formatter2 is None:
            tick_formatter2=FormatterDMS()

        klass = floating_axes.GridHelperCurveLinear
        klass.__init__(self, self._wcs_trans,
                       extremes,
                       grid_locator1=grid_locator1,
                       grid_locator2=grid_locator2,
                       tick_formatter1=tick_formatter1,
                       tick_formatter2=tick_formatter2)





class ParasiteAxesSky(ParasiteAxesAuxTrans):
    def __init__(self, parent, src_coord, *kl, **kwargs):

        dest_coord = parent.get_grid_helper()._wcsgrid_orig_coord_system
        tr = WcsSky2SkyTransform(src_coord, dest_coord) + \
             parent._wcsgrid_wcsaxes[0].transAux


        grid_helper = GridHelperWcsSky(parent.projection)
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
                          #tick_directiona="in",
                          axis_direction="bottom",
                          ):
        gh = self.get_grid_helper()
        axis = gh.new_floating_axis(nth_coord, value,
                                    #tick_direction=tick_direction,
                                    #label_direction=label_direction,
                                    axis_direction=axis_direction,
                                    axes=self._parent_axes)
        return axis


    def update_wcsgrid_params(self, **ka):
        """
        This method is deprecated. Use set_ticklabel_type or locator_params.
        """
        self.get_grid_helper().update_wcsgrid_params(**ka)



def get_transformed_image(Z, tr, extent=None, oversample=1.5):

    from scipy.ndimage import map_coordinates

    ny, nx = Z.shape
    if extent is None:
        #extent = 0, nx-1, 0, ny-1
        x1, x2, y1, y2 = 0, nx-1, 0, ny-1
    else:
        x1, x2, y1, y2 = extent

    transformed_edge = tr.transform(np.array([(x1, y1), (x1, y2), (x2, y1), (x2, y2)]))
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

    # this is a hack to find out the default interpolation mode for
    # imshow with affine. It was "nearest" but changed to "none".
    def get_imshow_affine_default_interpolation():
        from matplotlib.image import AxesImage

        class MyRenderer(object):
            def option_scale_image(self):
                return True

        class TestNearest(object):
            def get_interpolation(self):
                return "nearest"

        check_unsampled_image = AxesImage._check_unsampled_image.im_func
        if_nearest = check_unsampled_image(TestNearest(),
                                           MyRenderer())

        if if_nearest:
            return "nearest"
        else:
            return "none"

    _imshow_affine_default_interpolation = get_imshow_affine_default_interpolation()
    del get_imshow_affine_default_interpolation


    def imshow_affine(self, *kl, **kwargs):


        mode = self._imshow_affine_default_interpolation
        if kwargs.setdefault("interpolation", mode) != mode:
            raise ValueError("interpolation parameter must not be set or should be '%s'" % mode)
        im = ParasiteAxesAuxTrans.imshow(self, *kl, **kwargs)
        x1, x2, y1, y2 = im.get_extent()
        im._image_skew_coordinate = (x2, y1)

        return im



# def _parse_label_type(labtyp):
#     labtyp_parsed = labtyp.split(",")
#     labtyp_name = labtyp_parsed[0].strip()
#     labtyp_kwargs = {}
#     for p in labtyp_parsed[1:]:
#         kv = p.strip().split("=")
#         if len(kv) == 1:
#             labtyp_kwargs[kv]=True
#         else:
#             labtyp_kwargs[kv[0]]=kv[1]

#     return labtyp_name, labtyp_kwargs


# def _labtyp_updater_nbins(grid_helper, axis_num, v):
#     grid_finder = grid_helper.grid_finder
#     locator = [grid_finder.grid_locator1, grid_finder.grid_locator2][axis_num]
#     locator.den = int(v)

# _labtyp_updater_map = {}
# _labtyp_updater_default = dict(nbins=_labtyp_updater_nbins)

# def _labtyp_update(gh, axis_num, labtyp_name, k, v):
#     _labtyp_updater = _labtyp_updater_map.get(labtyp_name, _labtyp_updater_default)
#     _labtyp_updater[k](gh, axis_num, v)


# These methods need to go into the axisartist module.
class HostAxes(axisartist.HostAxes):
    def locator_params(self, axis='both', **kwargs):
        gh = self.get_grid_helper()
        gh.locator_params(axis=axis, **kwargs)
    locator_params.__doc__ = GridHelperCurveLinear.locator_params.__doc__


    def _set_label_text(self, axis, label, fontdict = None, **kwargs):
        axis.label.set_text(label)
        if fontdict is not None: axis.label.update(fontdict)
        axis.label.update(kwargs)
        return axis.label

    def set_xlabel(self, label, fontdict = None, **kwargs):
        self._set_label_text(self.axis["bottom"],
                             label, fontdict=fontdict, **kwargs)
        self._set_label_text(self.axis["top"],
                             label, fontdict=fontdict, **kwargs)

    def set_ylabel(self, label, fontdict = None, **kwargs):
        self._set_label_text(self.axis["left"],
                             label, fontdict=fontdict, **kwargs)
        self._set_label_text(self.axis["right"],
                             label, fontdict=fontdict, **kwargs)

    def get_xlabel(self, label):
        return self.axis["bottom"].label.get_text()

    def get_ylabel(self, label):
        return self.axis["left"].label.get_text()


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
                projection = get_kapteyn_projection(header)
            elif (header is None) and (wcs is not None):
                projection = get_kapteyn_projection(wcs)
            else:
                raise ValueError("wcs")

            self.projection = projection.sub([1,2])
            ctype1, ctype2 = self.projection.ctypes[:2]
            equinox = self.projection.equinox
            default_system = coord_system_guess(ctype1, ctype2, equinox)

            if default_system:
                grid_helper = GridHelperWcsSky(self.projection)
            else:
                grid_helper = GridHelperWcsSimple(self.projection)
            kw["grid_helper"] = grid_helper
        else:
            self.projection = kw["grid_helper"].projection

        super(AxesWcs, self).__init__(*kl, **kw)

        self._init_parasites()

        self.set_default_label("default", "default")

    def __iter__(self):
        raise TypeError("AxesWcs object is not iterable")


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
                proj = get_kapteyn_projection(key)
            except:
                raise ValueError("invalide key : %s" % repr(key))

            key = proj.sub([1,2])

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
        This method is deprecated. Use set_ticklabel_type or locator_params.

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

    def _get_default_label1(self, ctype_name, equinox):
        if ctype_name.startswith("GLON"):
            label=r"$l$"
        elif ctype_name.startswith("GLAT"):
            label=r"$b$"
        elif ctype_name.startswith("RA"):
            if equinox == 2000:
                label=r"$\alpha_{2000}$"
            elif equinox == 1950:
                label=r"$\alpha_{1950}$"
            else:
                label=r"$\alpha_{%6.1f}$" % (equinox)
        elif ctype_name.startswith("DEC"):
            if equinox == 2000:
                label=r"$\delta_{2000}$"
            elif equinox == 1950:
                label=r"$\delta_{1950}$"
            else:
                label=r"$\delta_{%6.1f}$" % (equinox)
        else:
            label = ""

        return label


    def _get_default_label_using_ctypes(self):
        ctype1_name, ctype2_name = self.projection.ctypes[:2]
        equinox = self.projection.equinox

        label1 = self._get_default_label1(ctype1_name, equinox)
        label2 = self._get_default_label1(ctype2_name, equinox)
        #coord_system = self.get_grid_helper().get_display_coord_system()

        return label1, label2
        #self.axis["right"].label.set_text(ylabel)
        #self.axis["bottom"].label.set_text(xlabel)
        #self.axis["top"].label.set_text(xlabel)
        #self.set_xlabel(xlabel)
        #self.set_ylabel(ylabel)
        #self.axis["bottom"].label.set_text(xlabel)
        #self.axis["top"].label.set_text(xlabel)

    def _get_default_label(self):
        coord_system = self.get_grid_helper().get_display_coord_system()
        if coord_system is None:
            coord_system = self.get_grid_helper()._wcsgrid_orig_coord_system

        if coord_system == "fk5":
            label1=r"$\alpha_{2000}$"
            label2=r"$\delta_{2000}$"
        elif coord_system == "fk4":
            label1=r"$\alpha_{1950}$"
            label2=r"$\delta_{1950}$"
        elif coord_system == "gal":
            label1=r"$l$"
            label2=r"$b$"
        else:
            label1=r""
            label2=r""

        return label1, label2


    def _decorate_default_label(self, label1, ticktyp1):
        if ticktyp1 in ["delta", "arcdeg", "arcmin", "arcsec", "arcmas"]:
            label1 = r"$\Delta$"+label1
            if ticktyp1 == "arcdeg":
                label1 = label1 + r" [$^{\circ}$]"
            elif ticktyp1 == "arcmin":
                label1 = label1 + r" [$^{\prime}$]"
            elif ticktyp1 == "arcsec":
                label1 = label1 + r" [$^{\prime\prime}$]"
            elif ticktyp1 == "arcmas":
                label1 = label1 + r" [mas]"

        else:
            if ticktyp1 == "absdeg":
                label1 = label1 + r" [$^{\circ}$]"
        return label1

    def set_default_label(self, ticktyp1=None, ticktyp2=None):
        """
        set x and y labels with the default values that are
        automatically determined according to the ticktyp and ctype of
        the axis.
        """

        if self.get_grid_helper()._wcsgrid_orig_coord_system is not None:
            label1, label2 = self._get_default_label()
        else:
            label1, label2 = self._get_default_label_using_ctypes()

        if ticktyp1 is not None:
            label1 = self._decorate_default_label(label1, ticktyp1)
            self.axis["bottom","top"].label.set_text(label1)

        if ticktyp2 is not None:
            label2 = self._decorate_default_label(label2, ticktyp2)
            self.axis["left","right"].label.set_text(label2)



    def set_ticklabel1_type(self, labtyp, **labtyp_kwargs):
        """
        Set the label type of the x axis. Available options
        depend on the grid_helper used. The default is
        GridHelperWcsBase and its variants.

        See
        :meth:`~pywcsgrid2.axes_wcs.GridHelperWcsBase.set_ticklabel1_type`
        for available options.
        """
        self.get_grid_helper().set_ticklabel1_type(labtyp, **labtyp_kwargs)

    def set_ticklabel2_type(self, labtyp, **labtyp_kwargs):
        """
        Set the label type of the y axis. Available options
        depend on the grid_helper used. The default is
        GridHelperWcsBase and its variants.

        See
        :meth:`~pywcsgrid2.axes_wcs.GridHelperWcsBase.set_ticklabel1_type`
        for available options.
        """
        self.get_grid_helper().set_ticklabel2_type(labtyp, **labtyp_kwargs)

    def set_ticklabel_type(self, labtyp1, labtyp2=None,
                           center_pixel=None,
                           labtyp1_kwargs=None,
                           labtyp2_kwargs=None,
                           **kwargs
                           ):
        """
        Set the label type of the x and y axes. Available options
        depend on the grid_helper used. The default is
        GridHelperWcsBase and its variants.

        See
        :meth:`~pywcsgrid2.axes_wcs.GridHelperWcsBase.set_ticklabel_type`
        for available options.
        """

        # Borrowed from Miriad's cgdisp.

        if labtyp2 is None:
            labtyp2 = labtyp1

        if center_pixel is None:
            # use the center of the current viewlim
            x, y, w, h = self.viewLim.bounds
            cx, cy = x + w/2., y + h/2.
        else:
            cx, cy = center_pixel

        if labtyp1_kwargs is None:
            labtyp1_kwargs = {}
        if labtyp2_kwargs is None:
            labtyp2_kwargs = {}
        labtyp1_kwargs.update(kwargs)
        labtyp2_kwargs.update(kwargs)

        gh = self.get_grid_helper()
        gh.set_ticklabel_type(labtyp1, labtyp2,
                              center_pixel=(cx, cy),
                              labtyp1_kwargs=labtyp1_kwargs,
                              labtyp2_kwargs=labtyp2_kwargs)
        self.set_default_label(labtyp1, labtyp2)


    def set_label_type(self, labtyp1, labtyp2,
                       offset_center=None):
        warnings.warn("set_label_type is deprecated. Use set_ticklabel_type")
        self.set_ticklabel_type_deprecated(labtyp1, labtyp2,
                                           offset_center=offset_center)

    # def set_label_type(self, labtyp1, labtyp2,
    #                    offset_center=None,
    #                    ):

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


    default_path_effects = []

    def set_default_path_effects(self, path_effects):
        """
        *default_path_effects* is applied to all the artists created by
        add-* methods. *default_path_effects* must be a list of path effects.
        """
        self.default_path_effects = path_effects

    def add_beam_size(self, major_pixel, minor_pixel, angle, loc,
                      frameon=False, borderpad=.8, patch_props=None,
                      path_effects=None,
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

        See mpl_toolkits.axes_grid1.anchored_artists for orther options.
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

        if path_effects is None:
            path_effects = self.default_path_effects

        if path_effects:
            ae.ellipse.set_path_effects(path_effects)

        self.add_artist(ae)
        return ae


    def add_inner_title(self, title, loc, path_effects=None, **kwargs):
        """
        Add a text at the inside corner of the axes.
        It simply uses mpl_toolkits.axes_grid1.anchored_artists.AnchoredText.
        See mpl_toolkits.axes_grid1.anchored_artists.AnchoredText
        for more details
        """
        # Figure title

        if path_effects is None:
            path_effects = self.default_path_effects

        if path_effects:
            kwargs.setdefault("frameon", False)
            at = AnchoredText(title, loc=loc, **kwargs)
            at.txt._text.set_path_effects(path_effects)
        else:
            at = AnchoredText(title, loc=loc, **kwargs)

        self.add_artist(at)

        return at


    def add_compass(self, loc, coord="fk5", arrow_length=0.15,
                    txt1="E", txt2="N",
                    delta_a1=0, delta_a2=0,
                    path_effects=None,
                    **kwargs):
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
        if path_effects is None:
            path_effects = self.default_path_effects

        if path_effects:
            ac.set_path_effects(path_effects)

        self.add_artist(ac)
        return ac


    def add_size_bar(self, length_pixel, label, loc, sep=5,
                     borderpad=0.8, frameon=False,
                     path_effects=None, **kwargs):
        """
        Add a horizontal line with a label underneath. Intended to
        display the angular size.

        Keyword arguments:

          *length_pixel* : length of the line in pixel in the image coordinate.

          *label* : label

          *loc* : the location code as in the legend

        """

        if path_effects is None:
            path_effects = self.default_path_effects

        asb =  AnchoredSizeBar(self.transData,
                               length_pixel,
                               label,
                               loc=loc,
                               borderpad=borderpad, sep=sep,
                               frameon=frameon,
                               **kwargs)
        if path_effects:
            for a in [asb.size_bar._children[0],
                      asb.txt_label._text]:
                a.set_path_effects(path_effects)

        self.add_artist(asb)

        return asb



SubplotWcs = maxes.subplot_class_factory(AxesWcs)






# class AxesWcsSky(AxesWcs):

#     def set_default_label2(self):
#         coord_system = self.get_grid_helper().get_display_coord_system()
#         if coord_system is None:
#             coord_system = self.get_grid_helper()._wcsgrid_orig_coord_system

#         if coord_system == "fk5":
#             label1=r"$\alpha_{2000}$"
#             label2=r"$\delta_{2000}$"
#         elif coord_system == "fk4":
#             label1=r"$\alpha_{1950}$"
#             label2=r"$\delta_{1950}$"
#         elif coord_system == "gal":
#             label1=r"$l$"
#             label2=r"$b$"

#         if label1 and self.get_grid_helper()._center_world is not None:
#             label1 = r"$\Delta$"+label1
#             label2 = r"$\Delta$"+label2

#         self.axis["bottom","top"].label.set_text(label1)
#         self.axis["left","right"].label.set_text(label2)







