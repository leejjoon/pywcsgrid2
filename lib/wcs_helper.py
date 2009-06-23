import pyfits
import kapteyn.wcs

class wcs(object):
    def __init__(self, header):
        if isinstance(header, pyfits.Header):
            self.projection = kapteyn.wcs.Projection(header)
        elif isinstance(header, kapteyn.wcs.Projection):
            self.projection = header
        elif hasattr(header, "wcs") and hasattr(header.wcs, "to_header"):
            # pywcs
            header = header.wcs.to_header()
            self.projection = kapteyn.wcs.Projection(header)
        else:
            raise ValueError("Unknown input type : %s", type(header))

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
