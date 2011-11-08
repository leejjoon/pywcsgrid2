import matplotlib.pyplot as plt
import pyfits

import pywcsgrid2 #.axes_wcs import GridHelperWcsFloating, AxesWcs


def demo_header():
    # header retrieved from "lambda_mollweide_halpha_fwhm06_0512.fits"
    header = """XTENSION= 'IMAGE   '           / IMAGE extension
BITPIX  =                  -32 / Number of bits per data pixel
NAXIS   =                    2 / Number of data axes
NAXIS1  =                  199 /
NAXIS2  =                   19 /
PCOUNT  =                    0 / No Group Parameters
GCOUNT  =                    1 / One Data Group
EXTNAME = 'TEMPERATURE'
CTYPE1  = 'RA---SIN'           / Coordinate Type
CTYPE2  = 'DEC--SIN'           / Coordinate Type
EQUINOX =              2000.00 / Equinox of Ref. Coord.
CDELT1  =     -0.0791293637247 / Degrees / Pixel
CDELT2  =      0.0791293637247 / Degrees / Pixel
CROTA2  =              0.00000 / Rotation Angle (Degrees)
CRPIX1  =                  100 / Reference Pixel in X
CRPIX2  =                   10 / Reference Pixel in Y
CRVAL1  =        67.0000000000 / Galactic longitude of reference pixel
CRVAL2  =         3.0000000000 / Galactic latitude of reference pixel
HISTORY PUTAST: Jun 17 14:36:35 2009 World Coordinate System parameters written
"""
    cards = pyfits.CardList()
    for l in header.split("\n"):
        #card = pyfits.Card()
        #card.fromstring(l.strip())
        card = pyfits.Card.fromstring(l.strip())
        cards.append(card)
    h = pyfits.Header(cards)
    return h

h = demo_header()
nx, ny = h["naxis1"], h["naxis2"],

labtypes = [("hms", {}),
            ("h", {}),
            ("dms", {}),
            ("absdeg", {}),
            ("delta", dict(offset=h["crval1"], latitude=h["crval2"])),
            ("arcmin", dict(offset=h["crval1"], latitude=h["crval2"])),
            ("manual", dict(locs=[70, 65, 60], labels=["A","B","C"])),
            ]
n = len(labtypes)

plt.figure(figsize=(5, 1*n))

for i, (labtyp1, labtyp1_kwargs) in enumerate(labtypes):
    ax = pywcsgrid2.subplot(n, 1, i+1, header=h, frameon=False)
    ax.set_aspect(1.)
    ax.set_xlim(-0.5, nx-0.5)
    ax.set_ylim(-0.5, 4.5)
    ax.axis["left","right","top"].set_visible(False)

    ax.annotate(labtyp1, (0.,0), (0, 2),
                xycoords="axes fraction",
                textcoords="offset points")
    ax.set_ticklabel1_type(labtyp1, **labtyp1_kwargs)
    ax.set_default_label(labtyp1, None)
    

plt.show()
