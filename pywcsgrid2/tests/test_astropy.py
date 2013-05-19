from __future__ import absolute_import
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pywcsgrid2
#from cStringIO import StringIO
import io

def get_header():
    from pywcsgrid2.allsky_axes import allsky_header
    hdr = allsky_header("gal", "CAR", cdelt=1)
    return repr(hdr.ascard)

header_str = get_header()

def _demo_astropy(pyfits, pywcs=None):

    header = pyfits.Header.fromstring(header_str)

    plt.figure()
    if pywcs is None:
        ax = pywcsgrid2.subplot(111, header=header)
    else:
        wcs = pywcs.WCS(header)
        ax = pywcsgrid2.subplot(111, wcs=wcs)

    ax.imshow([[0, 0], [0,0]], origin="low", interpolation="none")
    ax.grid(color="w")

    s = io.BytesIO()
    plt.savefig(s, format="png")


def test_astropy():
    import astropy.io.fits
    import astropy.wcs
    _demo_astropy(astropy.io.fits)
    _demo_astropy(astropy.io.fits, astropy.wcs)

def test_pyfits_pywcs():
    import pyfits
    _demo_astropy(pyfits)
    try:
        import pywcs
    except ImportError:
        pass
    else:
        _demo_astropy(pyfits, pywcs)

if __name__ == "__main__":
    test_astropy()
    test_pyfits_pywcs()
