import matplotlib
from matplotlib.pyplot import gcf, gca, gci, draw_if_interactive

## Plotting part 1: manually generated functions and wrappers ##

import pywcsgrid2.axes_grid as axes_grid

def colorbar(mappable=None, cax=None, ax=None, **kw):
    if mappable is None:
        mappable = gci()

    if ax is None:
        ax = gca()

    ret = axes_grid.colorbar.colorbar(mappable, cax=cax, ax=ax, **kw)
    draw_if_interactive()
    return ret

colorbar.__doc__ = matplotlib.colorbar.colorbar_doc

