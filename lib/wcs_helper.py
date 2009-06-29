#import pyfits
#import kapteyn.wcs
import numpy as np

# class wcs(object):
#     def __init__(self, header):
#         if isinstance(header, pyfits.Header):
#             self.projection = kapteyn.wcs.Projection(header)
#         elif isinstance(header, kapteyn.wcs.Projection):
#             self.projection = header
#         elif hasattr(header, "wcs") and hasattr(header.wcs, "to_header"):
#             # pywcs
#             header = header.wcs.to_header()
#             self.projection = kapteyn.wcs.Projection(header)
#         else:
#             raise ValueError("Unknown input type : %s", type(header))

# def meshgrid_pcolor(wcs, img_shape):
#     ny, nx = img_shape

#     gx = np.arange(-0.5, nx+1., 1.)
#     gy = np.arange(-0.5, ny+1., 1.)

#     mx, my = meshgrid(gx, gy)

#     wx, wy = wcs.pixel2world(mx.flat, my.flat)
#     return wx.reshape([nx+1, ny+1]), wy.reshape([nx+1, ny+1])


# def meshgrid_contour(wcs, img_shape):
#     ny, nx = img_shape

#     gx = np.arange(0., nx+.5, 1.)
#     gy = np.arange(0., ny+.5, 1.)

#     mx, my = meshgrid(gx, gy)

#     wx, wy = wcs.pixel2world(mx.flat, my.flat)
#     return wx.reshape([nx, ny]), wy.reshape([nx, ny])


# def coord_system_guess(ctype1_name, ctype2_name, equinox):
#     return None

# class coord_sys(object):
#     pass





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

