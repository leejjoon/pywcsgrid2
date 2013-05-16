import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pywcsgrid2
from cStringIO import StringIO

def _demo_astropy(pyfits, pywcs=None):
    f = pyfits.open("data/lmc.fits")
    plt.figure()
    if pywcs is None:
        ax = pywcsgrid2.subplot(111, header=f[0].header)
    else:
        wcs = pywcs.WCS(f[0].header)
        ax = pywcsgrid2.subplot(111, wcs=wcs)

    ax.imshow(f[0].data, origin="low",
              vmin=0, vmax=2000, interpolation="none")
    ax.grid(color="w")

    s = StringIO()
    plt.savefig(s, format="png")


def test_astropy():
    import astropy.io.fits
    import astropy.wcs
    _demo_astropy(astropy.io.fits)
    _demo_astropy(astropy.io.fits, astropy.wcs)

    import pyfits
    import pywcs
    _demo_astropy(pyfits)
    _demo_astropy(pyfits, pywcs)
