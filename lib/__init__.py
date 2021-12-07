from __future__ import absolute_import

from .axes_wcs import SubplotWcs as Subplot
from .axes_wcs import AxesWcs as Axes
#from axes_wcs import GridHelperWcs as GridHelper
from .axes_wcs import GridHelperWcsSky, GridHelperWcsSimple

GridHelperSky, GridHelperSimple = GridHelperWcsSky, GridHelperWcsSimple
GridHelper = GridHelperSky

from .floating_axes import FloatingAxes, FloatingSubplot

import matplotlib.pyplot as plt

def _check_kwargs(kwargs):
    if ("header" in kwargs) or ("grid_helper" in kwargs) or ("wcs" in kwargs):
        return

    raise ValueError("one of header, wcs, or grid_helper parameter must be specified " )


def axes(*args, **kwargs):
    """
    Creates an axes for displaying FITS images. The axes uses the
    FITS header information for ticks and grids appropriate for
    sky coordinate.

    The arguments for this axes is same as mpl's Axes, except that
    either *header* or *grid_helper* (but not both) keyword
    arguments must be provided.

    Necessary Keyword arguments:

    *grid_helper*: pywcsgrid2.GridHelper instance

    or:

    *header*: pyfits.Header instance header of the fits file that will
              be used for the axes. A new pywcsgrid2.GridHelper
              instance is created using the header.

    *axis_nums*: a tuple of two integers specifying which axis of fits file
                 will be used for x- and y-axis of image. Default is (0, 1).
                 If CTYPE1=="DEC" and CTYPE2=="RA", you may use (1, 0).
                 This can also be used for n-d image.

    """
    _check_kwargs(kwargs)

    nargs = len(args)
    if len(args)==0: return subplot(111, **kwargs)
    if nargs>1:
        raise TypeError('Only one non keyword arg to axes allowed')
    arg = args[0]

    if isinstance(arg, Axes):
        a = plt.gcf().sca(arg)
    else:
        rect = arg
        ax = Axes(plt.gcf(), rect, **kwargs)
        a = plt.gcf().add_axes(ax)
    plt.draw_if_interactive()
    return a


def subplot(*args, **kwargs):
    """
    """

    _check_kwargs(kwargs)

    fig = plt.gcf()
    ax = Subplot(fig, *args, **kwargs)
    fig.add_subplot(ax)
    return ax


def floating_axes(*args, **kwargs):

    assert("header" in kwargs)
    assert("extremes" in kwargs)

    nargs = len(args)
    if len(args)==0: return floating_subplot(111, **kwargs)
    if nargs>1:
        raise TypeError('Only one non keyword arg to axes allowed')
    arg = args[0]

    if isinstance(arg, Axes):
        a = plt.gcf().sca(arg)
    else:
        rect = arg
        ax = FloatingAxes(plt.gcf(), rect, **kwargs)
        a = plt.gcf().add_axes(ax)
    plt.draw_if_interactive()
    return a


def floating_subplot(*args, **kwargs):
    """
    """

    assert("header" in kwargs)
    assert("extremes" in kwargs)

    fig = plt.gcf()
    ax = FloatingSubplot(fig, *args, **kwargs)
    ax = plt.subplot(ax)
    return ax


import re

_wcs_key_pattern = re.compile(r'^(NAXIS|CD|CDELT|CRPIX|CRVAL|CTYPE|CROTA|LONGPOLE|LATPOLE|PV|DISTORT|OBJECT|BUNIT|EPOCH|EQUINOX|LTV|LTM|DTV|DTM)')

def filterwcs(h):
    """
    select wcs related cards
    """
    from astropy_helper import pyfits
    # We have to re-instantiate the cards since we don't know if the original
    # header was from PyFITS or Astropy.
    l = [pyfits.Card(card.key, card.value, card.comment) for card in h.ascardlist() if _wcs_key_pattern.match(card.key)]
    return pyfits.Header(l)


