import numpy as np

from .astropy_helper import pyfits, pywcs

import healpy
import warnings


class HealpixData(object):
    def __init__(self, nside, data, coord=None, nested=False, flipy=False):
        self._nside = nside
        self._data = data
        self._nested = nested
        self._flipy = flipy
        self._coord = coord

    def get_projected_map(self, header):

        map_shape = (header["naxis2"], header["naxis1"])
        iy, ix = np.indices(map_shape)
        wcs = pywcs.WCS(header)
        phi, theta = wcs.wcs_pix2world(ix, iy, 0)

        if self._coord is not None:
            from pywcsgrid2.wcs_helper import coord_system_guess, sky2sky
            map_coord = coord_system_guess(header["ctype1"],
                                           header["ctype2"],
                                           equinox=header["equinox"])
            if (map_coord is not None) and (map_coord != self._coord):
                warnings.warn(" doing the conversion " + map_coord)
                phi, theta = sky2sky(map_coord, self._coord)(phi, theta)


        if self._flipy:
            theta -= 90
            theta *= -np.pi/180.
        else:
            theta += 90
            theta *= np.pi/180.

        phi *= np.pi/180

        if self._nested:
            ang2pix = healpy._healpy_pixel_lib._ang2pix_nest
        else:
            ang2pix = healpy._healpy_pixel_lib._ang2pix_ring


        # some values could be NaNs. Maske those out before calling
        # ang2pix and recover them.
        mask = np.isfinite(theta) & np.isfinite(theta)

        ipix = ang2pix(self._nside, theta[mask], phi[mask])

        map_data_ = self._data[ipix]
        map_data = np.empty(map_shape, dtype=map_data_.dtype)
        map_data.fill(np.nan)
        map_data[mask] = map_data_

        return map_data


if __name__ == '__main__':

    fname = "LAB_fullvel.fits"
    f = pyfits.open(fname)
    header = f[1].header

    ordering = header["ordering"]
    nside = header["nside"]
    data = f[1].data["temperature"]

    healpix_data = HealpixData(nside, data.flat, nested=False)

    fits_name = "lambda_mollweide_halpha_fwhm06_0512.fits"
    f2 = pyfits.open(fits_name)

    d = healpix_data.get_projected_map(f2[1].header)

    #data2 = f2[1].data
    #header2 = f2[1].header



