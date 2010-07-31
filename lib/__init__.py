
from axes_wcs import SubplotWcs as Subplot
from axes_wcs import AxesWcs as Axes
#from axes_wcs import GridHelperWcs as GridHelper
from axes_wcs import GridHelperWcsSky, GridHelperWcsSimple

GridHelperSky, GridHelperSimple = GridHelperWcsSky, GridHelperWcsSimple
GridHelper = GridHelperSky

import matplotlib.pyplot as plt

def _check_kwargs(kwargs):
    if ("header" in kwargs) or ("grid_helper" in kwargs) or ("wcs" in kwargs):
        return

    raise ValueError("one of header, wcs, or grid_helper parameter must be specified " )


def axes(*args, **kwargs):

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
    ax = plt.subplot(ax)
    return ax


import re

_wcs_key_pattern = re.compile(r'^(NAXIS|CD|CDELT|CRPIX|CRVAL|CTYPE|CROTA|LONGPOLE|LATPOLE|PV|DISTORT|OBJECT|BUNIT|EPOCH|EQUINOX|LTV|LTM|DTV|DTM)')

def filterwcs(h):
    """
    select wcs related cards
    """
    import pyfits
    l = [card for card in h.ascardlist() if _wcs_key_pattern.match(card.key)]
    return pyfits.Header(l)



