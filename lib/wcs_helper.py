import numpy as np

from kapteyn_celestial import skymatrix, longlat2xyz, dotrans, xyz2longlat
import kapteyn_celestial

import pyfits

_wcs_module_import_log = []

_pywcs_installed = False
try:
    import pywcs
except ImportError:
    _wcs_module_import_log.append("Failed to import the pywcs")
else:
    if hasattr(pywcs.WCS, "sub"):
        _pywcs_installed = True
    else:
        _wcs_module_import_log.append("pywcs imported but does not have 'sub' attribute. More recent version of pywcs is required.")

_kapteyn_installed = False
try:
    import kapteyn.wcs
except ImportError:
    _wcs_module_import_log.append("Failed to import the kpateyn.wcs")
else:
    _kapteyn_installed = True

if not _kapteyn_installed and not _pywcs_installed:
    err = ["Either pywcs or Kapteyn python packages are required."]
    err.extend(_wcs_module_import_log)

    raise ImportError("\n".join(err))

FK4 = (kapteyn_celestial.equatorial, kapteyn_celestial.fk4, 'B1950.0')
FK5 = (kapteyn_celestial.equatorial, kapteyn_celestial.fk5, 'J2000.0')
GAL = kapteyn_celestial.galactic

coord_system = dict(fk4=FK4,
                    fk5=FK5,
                    gal=GAL)

def is_equal_coord_sys(src, dest):
    return (src.lower() == dest.lower())

def is_string_like(obj):
    'Return True if *obj* looks like a string'
    if isinstance(obj, (str, unicode)): return True
    try: obj + ''
    except: return False
    return True

class sky2sky(object):
    def __init__(self, src, dest):

        if is_string_like(src):
            src = coord_system[src.lower()]

        if is_string_like(dest):
            dest = coord_system[dest.lower()]

        self.src = src
        self.dest = dest
        self._skymatrix = skymatrix(src, dest)
        #self.tran = wcs.Transformation(src, dest)

    def inverted(self):
        return sky2sky(self.dest, self.src)

    def _dotran(self, lonlat):
        xyz = longlat2xyz(lonlat)
        xyz2 = dotrans(self._skymatrix, xyz)
        lonlats2 = xyz2longlat(xyz2)
        return lonlats2

    def __call__(self, lon, lat):
        lon, lat = np.asarray(lon), np.asarray(lat)
        lonlat = np.concatenate([lon[:,np.newaxis],
                                 lat[:,np.newaxis]],1)
        ll_dest = np.asarray(self._dotran(lonlat))
        return ll_dest[:,0], ll_dest[:,1]


def coord_system_guess(ctype1_name, ctype2_name, equinox):
    if ctype1_name.upper().startswith("RA") and \
       ctype2_name.upper().startswith("DEC"):
        if equinox == 2000.0:
            return "fk5"
        elif equinox == 1950.0:
            return "fk4"
        elif equinox is None:
            return "fk5"
        else:
            return None
    if ctype1_name.upper().startswith("GLON") and \
       ctype2_name.upper().startswith("GLAT"):
        return "gal"
    return None


class ProjectionBase(object):
    """
    A wrapper for kapteyn.projection or pywcs
    """
    def _get_ctypes(self):
        pass

    def _get_equinox(self):
        pass

    def topixel(self):
        """ 1, 1 base """
        pass

    def toworld(self):
        """ 1, 1 base """
        pass

    def sub(self, axes):
        pass


class ProjectionKapteyn(ProjectionBase):
    """
    A wrapper for kapteyn.projection 
    """
    def __init__(self, header):
        if isinstance(header, pyfits.Header):
            self._proj = kapteyn.wcs.Projection(header)
        else:
            self._proj = header

    def _get_ctypes(self):
        return self._proj.ctype

    ctypes = property(_get_ctypes)

    def _get_equinox(self):
        return self._proj.equinox

    equinox = property(_get_equinox)

    def topixel(self, xy):
        """ 1, 1 base """
        return self._proj.topixel(xy)

    def toworld(self, xy):
        """ 1, 1 base """
        return self._proj.toworld(xy)

    def sub(self, axes):
        proj = self._proj.sub(axes=axes)
        return ProjectionKapteyn(proj)


class ProjectionPywcs(ProjectionBase):
    """
    A wrapper for pywcs
    """
    def __init__(self, header):
        if isinstance(header, pyfits.Header):
            self._pywcs = pywcs.WCS(header=header)
        else:
            self._pywcs = header

    def _get_ctypes(self):
        return tuple(self._pywcs.wcs.ctype)

    ctypes = property(_get_ctypes)

    def _get_equinox(self):
        return self._pywcs.wcs.equinox

    equinox = property(_get_equinox)

    def topixel(self, xy):
        """ 1, 1 base """
        xy2 = self._pywcs.wcs_sky2pix(np.asarray(xy).T, 1)
        return xy2[:,0], xy2[:,1]

    def toworld(self, xy):
        """ 1, 1 base """
        xy2 = self._pywcs.wcs_pix2sky(np.asarray(xy).T, 1)
        return xy2[:,0], xy2[:,1]

    def sub(self, axes):
        wcs = self._pywcs.sub(axes=axes)
        return ProjectionPywcs(wcs)


if _pywcs_installed:
    ProjectionDefault = ProjectionPywcs
else:
    ProjectionDefault = ProjectionKapteyn

def get_kapteyn_projection(header):
    if _kapteyn_installed and isinstance(header, kapteyn.wcs.Projection):
        projection = ProjectionKapteyn(header)
    elif _pywcs_installed and isinstance(header, pywcs.WCS):
        projection = ProjectionPywcs(header)
    elif isinstance(header, ProjectionBase):
        projection = header
    else:
        projection = ProjectionDefault(header)

    projection = projection.sub(axes=[1,2])
    return projection


def estimate_cdelt(transSky2Pix, x0, y0):

    transPix2Sky = transSky2Pix.inverted()

    lon0, lat0 = transPix2Sky.transform_point((x0, y0))

    lon1, lat1 = transPix2Sky.transform_point((x0+1, y0))
    dlon = (lon1-lon0)*np.cos(lat0/180.*np.pi)
    dlat = (lat1-lat0)
    cd1 = (dlon**2 + dlat**2)**.5

    lon2, lat2 = transPix2Sky.transform_point((x0+1, y0))
    dlon = (lon2-lon0)*np.cos(lat0/180.*np.pi)
    dlat = (lat2-lat0)
    cd2 = (dlon**2 + dlat**2)**.5

    return (cd1*cd2)**.5


def estimate_angle(transSky2Pix, x0, y0):
    """
    return a tuple of two angles (in degree) of increasing direction
    of 1st and 2nd coordinates.

    note that x, y = wcs_proj.topixel(sky_to_sky((l1, l2)))

    """

    cdelt = estimate_cdelt(transSky2Pix, x0, y0)

    transPix2Sky = transSky2Pix.inverted()

    lon0, lat0 = transPix2Sky.transform_point((x0, y0))

    x1, y1 = transSky2Pix.transform_point((lon0 + cdelt*np.cos(lat0/180.*np.pi),
                                           lat0))

    x2, y2 = transSky2Pix.transform_point((lon0, lat0+cdelt))

    a1 = np.arctan2(y1-y0, x1-x0)/np.pi*180.
    a2 = np.arctan2(y2-y0, x2-x0)/np.pi*180.

    return a1, a2




if __name__ == "__main__":
    fk5_to_fk4 = sky2sky(FK5, FK4)
    print fk5_to_fk4([47.37], [6.32])
    print fk5_to_fk4([47.37, 47.37], [6.32, 6.32])
    print sky2sky("fk5", "FK4")([47.37, 47.37], [6.32, 6.32])




