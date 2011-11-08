import matplotlib.pyplot as plt
import pyfits

from mpl_toolkits.axisartist.floating_axes import floatingaxes_class_factory

from pywcsgrid2.axes_wcs import GridHelperWcsFloating, AxesWcs

import matplotlib.axes as maxes


# FIXME : It would be possible to directly create WCS instance, instead of
# using pyfit sheader. -JJL


# these sizes are for cdelt = 0.1. naxis needs to be twice of these values
# For MER projection, the y size is for 0-75 degree.

_allsky_image_size_map = dict(PAR=(1837, 919),
                              CYP=(1837, 1170),
                              CEA=(1837, 586),
                              MOL=(1654, 828),
                              CAR=(1837, 919),
                              MER=(1837, 1186),
                              SFL=(1837, 919))

def allsky_header(coord, proj, lon_center=0., cdelt=0.1):
    # header retrieved from "lambda_mollweide_halpha_fwhm06_0512.fits"

    proj = proj.upper()
    nx, ny =_allsky_image_size_map.get(proj, (4096, 2048))

    from math import ceil
    nx, ny = int(ceil(nx/cdelt*0.1)), int(ceil(ny/cdelt*0.1))

    header = """XTENSION= 'IMAGE   '           / IMAGE extension
BITPIX  =                  -32 / Number of bits per data pixel
NAXIS   =                    2 / Number of data axes
NAXIS1  =                 %d /
NAXIS2  =                 %d /
CDELT1  =     -%f / Degrees / Pixel
CDELT2  =      %f / Degrees / Pixel
CROTA2  =              0.00000 / Rotation Angle (Degrees)
CRPIX1  =              %.1f / Reference Pixel in X
CRPIX2  =              %.1f / Reference Pixel in Y
CRVAL2  =        0.00000000000 / Galactic latitude of reference pixel
LONPOLE =      180.00000000000 / Native longitude of Galactic pole
LATPOLE =       90.00000000000 / Galactic latitude of native pole
""" % (2*nx, 2*ny, cdelt, cdelt, nx+0.5, ny+0.5)

    header_list = header.split("\n")
    headeradd = header_list.append
    
    if coord.startswith("fk"):
        headeradd("CTYPE1  = 'RA---%s'           / Coordinate Type" % proj)
        headeradd("CTYPE2  = 'DEC--%s'           / Coordinate Type" % proj)
        if coord == "fk5":
            headeradd("EQUINOX =              2000.00 / Equinox of Ref. Coord.")
        elif coord == "fk4":
            headeradd("EQUINOX =              1950.00 / Equinox of Ref. Coord.")
        else:
            raise ValueError()
    elif coord == "gal":
        headeradd("CTYPE1  = 'GLON-%s'           / Coordinate Type" % proj)
        headeradd("CTYPE2  = 'GLAT-%s'           / Coordinate Type" % proj)
        headeradd("EQUINOX =              2000.00 / Equinox of Ref. Coord.")
    else:
        raise ValueError()


    headeradd("CRVAL1  =  %10.5f / Galactic longitude of reference pixel" \
              % (lon_center,))
    
    cards = pyfits.CardList()
    for l in header_list:
        # check if fromstring is a classmethod
        if type(pyfits.Card.__dict__["fromstring"]) is classmethod:
            card = pyfits.Card.fromstring(l.strip())
        else:
            card = pyfits.Card()
            card.fromstring(l.strip())
        cards.append(card)

    h = pyfits.Header(cards)
    return h


FloatingAxes = floatingaxes_class_factory(AxesWcs)
FloatingSubplot = maxes.subplot_class_factory(FloatingAxes)


_proj_pseudo_cyl_list = ["SFL", "PAR", "MOL", "AIT"]
_proj_lat_limits = dict(MER= 75)

def make_allsky_axes_from_header(fig, rect, header, lon_center,
                                 lat_minmax=None, pseudo_cyl=None):
    

    proj = header["CTYPE1"].split("-")[-1]
    if pseudo_cyl is None:
        if proj in _proj_pseudo_cyl_list:
            pseudo_cyl = True
        
    if lat_minmax is None:
        lat_max = _proj_lat_limits.get(proj, 90)
        lat_min = -lat_max
    else:
        lat_min, lat_max = lat_minmax

    extremes = (lon_center+180., lon_center-180., lat_min, lat_max)
    
    grid_helper = GridHelperWcsFloating(wcs=header, extremes=extremes)

    ax = FloatingSubplot(fig, rect, grid_helper=grid_helper)
    ax.locator_params(axis="x", nbins=7)
    ax.locator_params(axis="y", nbins=8)
    ax.set_autoscale_on(False)
    grid_helper.set_lon_ref(lon_center-180)
    fig.add_subplot(ax)

    if pseudo_cyl:
        ax.axis["bottom", "top"].set_visible(False)

        axis = ax.axis["left"]
        axis.major_ticklabels._text_follow_ref_angle=False
        axis.major_ticklabels.set(rotation=0,
                              va="center", ha="center",
                              pad=0)

        ax.axis["lat=0"] = grid_helper.new_floating_axis(nth_coord=1, value=0,
                                                       axes=ax,
                                                       axis_direction='bottom')

        axis = ax.axis["lat=0"]
        axis.get_helper().set_extremes(lon_center+180., lon_center-180)
        axis.set_ticklabel_direction("-")
        axis.set_axislabel_direction("-")
        axis.label.set_text(ax.axis["bottom"].label.get_text())

    else:
        ax.axis["top"].set_ticklabel_direction("+")
        ax.axis["top"].set_axislabel_direction("+")
        ax.axis["bottom"].set_ticklabel_direction("-")
        ax.axis["bottom"].set_axislabel_direction("-")

        ax.axis["bottom"].major_ticklabels.set_pad(6)


    return ax



def make_allsky_axes(fig, rect, coord, proj, lon_center=0,
                     lat_minmax=None, pseudo_cyl=None):

    header = allsky_header(coord, proj, lon_center)

    ax = make_allsky_axes_from_header(fig, rect, header,
                                      lon_center,
                                      pseudo_cyl=pseudo_cyl)

    return ax


if __name__ == '__main__':
    
    proj_list = ["CYP", "CEA", "CAR", "MER", "SFL", "PAR", "MOL", "AIT" ]


    for proj in proj_list:
        fig = plt.figure()
        rect = 111
        ax = make_allsky_axes(fig, rect, "gal", proj, lon_center=0)
        ax.set_title("proj = %s" % proj, position=(0.5, 1.1))
    
        plt.show()
