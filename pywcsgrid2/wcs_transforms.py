from matplotlib.transforms import Transform
from matplotlib.path import Path
import numpy as np

# pixel origin need to be corrected!

from wcs_helper import coord_system_guess, sky2sky, get_kapteyn_projection

def coord_conv(src, dest):
    return sky2sky(src, dest)

class CurvedTransform(Transform):
    def __init__(self, resolution=None):
        """
        Create a new WCS transform.
        """
        Transform.__init__(self)
        if resolution is None:
            resolution = 1
        self._resolution = resolution

    def transform_path(self, path):
        vertices = path.vertices
        #ipath = path.interpolated(self._resolution)
        #return Path(self.transform(ipath.vertices), ipath.codes)
        return Path(self.transform(path.vertices), path.codes)

    transform_path.__doc__ = Transform.transform_path.__doc__

    transform_path_non_affine = transform_path
    transform_path_non_affine.__doc__ = Transform.transform_path_non_affine.__doc__


class WcsSky2PixelTransform(CurvedTransform):
    """
    """
    input_dims = 2
    output_dims = 2
    is_separable = False

    def __init__(self, header_or_wcs, src_coord=None, resolution=None):
        """
        Create a new WCS transform.
        """
        CurvedTransform.__init__(self, resolution)

        self.projection = get_kapteyn_projection(header_or_wcs)

        if src_coord is not None:
            self.src_coord = src_coord

            ctype1, ctype2 = self.projection.ctypes[:2]
            equinox = self.projection.equinox
            coord_guess = coord_system_guess(ctype1, ctype2, equinox)

            self.dest_coord = coord_guess
            if src_coord == coord_guess:
                self.coord_conv_func = None
            else:
                self.coord_conv_func = coord_conv(src_coord, coord_guess)
        else:
            self.src_coord = None
            self.dest_coord = None
            self.coord_conv_func = None


    def transform(self, ll):

        lon, lat = ll[:,0], ll[:,1]

        if self.coord_conv_func is not None:
            lon, lat = self.coord_conv_func(lon, lat)

        x, y = self.projection.topixel((lon, lat))
        x, y = x-1, y-1 # => origin=0,0
        xy = np.concatenate((x[:,np.newaxis], y[:,np.newaxis]), 1)

        return xy
        #origin=0

        #return self.wcs.wcs.s2p(ll, origin)["pixcrd"]
    #return self.wcs.wcs_sky2pix(ll, =1)


    transform.__doc__ = Transform.transform.__doc__

    transform_non_affine = transform
    transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__


    def inverted(self):
        return WcsPixel2SkyTransform(self.projection, #wcs,
                                     dest_coord=self.src_coord,
                                     resolution=self._resolution)
    inverted.__doc__ = Transform.inverted.__doc__


class WcsPixel2SkyTransform(CurvedTransform):
    """
    The base Aitoff transform.
    """
    input_dims = 2
    output_dims = 2
    is_separable = False


    def __init__(self, header_or_wcs, dest_coord=None, resolution=None):
        """
        Create a new WCS transform.  Resolution is the number of steps
        to interpolate between each input line segment to approximate its
        path.
        """
        CurvedTransform.__init__(self, resolution)

        #if wcs.naxis > 2:
        #    wcs = wcs.wcssub([1, 2])

        self.projection = get_kapteyn_projection(header_or_wcs)
        #self.wcs = wcs

        if dest_coord is not None:
            ctype1, ctype2 = self.projection.ctypes[:2]
            equinox = self.projection.equinox
            coord_guess = coord_system_guess(ctype1, ctype2, equinox)

            self.src_coord = coord_guess
            self.dest_coord = dest_coord
            if dest_coord == coord_guess:
                self.coord_conv_func = None
            else:
                self.coord_conv_func = coord_conv(coord_guess, dest_coord)
        else:
            self.src_coord = None
            self.dest_coord = None
            self.coord_conv_func = None



    def transform(self, xy):
        xy = np.asarray(xy)
        x, y = xy[:,0], xy[:,1]

        x, y = x+1, y+1 # => origin=1,1
        lon, lat = self.projection.toworld((x, y))

        if self.coord_conv_func is not None:
            lon, lat = self.coord_conv_func(lon, lat)

        ll = np.concatenate((lon[:,np.newaxis], lat[:,np.newaxis]), 1)

        return ll



    transform.__doc__ = Transform.transform.__doc__

    transform_non_affine = transform
    transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__


    def inverted(self):
        return WcsSky2PixelTransform(self.projection, #self.wcs,
                                     src_coord=self.dest_coord,
                                     resolution=self._resolution)
    inverted.__doc__ = Transform.inverted.__doc__



class WcsSky2SkyTransform(CurvedTransform):
    """
    The base Aitoff transform.
    """
    input_dims = 2
    output_dims = 2
    is_separable = False

    def __init__(self, src_coord, dest_coord, resolution=None):
        """
        Create a new WCS transform.  Resolution is the number of steps
        to interpolate between each input line segment to approximate its
        path.
        """

        CurvedTransform.__init__(self, resolution)

        self.src_coord = src_coord
        self.dest_coord = dest_coord

        self.coord_conv_func = coord_conv(src_coord, dest_coord)


    def transform(self, ll):
        lon1, lat1 = ll[:,0], ll[:,1]
        lon2, lat2 = self.coord_conv_func(lon1, lat1)

        return np.concatenate((lon2[:,np.newaxis], lat2[:,np.newaxis]),
                              1)


    transform.__doc__ = Transform.transform.__doc__

    transform_non_affine = transform
    transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__


    def inverted(self):
        return WcsSky2SkyTransform(self.dest_coord, self.src_coord)
    inverted.__doc__ = Transform.inverted.__doc__


