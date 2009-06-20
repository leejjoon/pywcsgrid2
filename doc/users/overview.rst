========
Overview
========

The pywcsgrid2 is a collection of helper classes to display
astronomical fits images in matplotlib. Its main functionality is to
draw ticks, ticklabels, and grids in an appropriate sky
coordinates. It is based on ``pgsbox`` routine included in the
``wcslib``.  However, unlike the the original ``pgsbox``, it does not
depends on the ``pgplot`` library .

.. contents::
   :depth: 1
   :local:


For this module to be useful, you need pyfits and pywcs installed.
You read in the fits file using pyfits and create a pywcs.WCS instance
using its header information. Then you create a subplot (or axes) to
display the data with the wcs information. ::

    import pyfits
    import pywcs
    import matplotlib.pyplot as plt
    from pywcsgrid2.axes_wcs import AxesWcs

    f = pyfits.open("data/lmc.fits")
    wcs = pywcs.WCS(f[0].header)

    fig = plt.figure(1, [5,5])
    ax = AxesWcs(fig, [0.2, 0.15, 0.7, 0.8], wcs=wcs)
    fig.add_axes(ax)

The data coordinate of the AxesWcs axes is the image (pixel) coordinate
((0-based!) of the fits file. But The ticklabels are displayed in the
sky coordinates. For a simple image display, you can do ::

    ax.imshow(f[0].data)


.. plot:: figures/demo_basic1.py



Again, the data coordinate of the AxesWcs axes is the image coordinate.
For example, xlim and ylim needs to be in image coordinates (0-based). ::

    # viewlimits in image coordinate
    ax.set_xlim(6, 81)
    ax.set_ylim(23, 98)


As in the MPL, the ``grid`` method draws grid lines. It will draw
curved grid lines using the given wcs information. Note that
the ticks are also rotated accordingly.::

    ax.grid()

You can not use the mpl methods like ``set_ticks``. Instead, you need
to change the wcsgrid parameters associated. For example, to set the
approximate number of ticks in each axis, ::

  # change grid density
  gh = ax.get_grid_helper()
  gh.update_wcsgrid_params(label_density=[4,3])

.. plot:: figures/demo_basic2.py

Or explicitly set cooridnates (in degrees) for ticks and grids.::

  gh.update_wcsgrid_params(grid2=[-65, -70, -75])

You can change the displyed sky coordinate (i.e., coordinates for
ticks, ticklabels and grids). For example, to display the Galactic
coordinate system::

    ax.set_display_coord_system("gal")

Sometimes, you will need to swap the axis for better tick labeling
(i.e., xaxis display latitude and yaxis display longitude). ::

    ax.swap_tick_coord()

The AxesWcs class is based on the axes_grid.axislines.Axes. For
eaxmple, to turn on the top and right tick labels,::

  ax.axis["top"].major_ticklabels.set_visible(True)
  ax.axis["right"].major_ticklabels.set_visible(True)

.. plot:: figures/demo_basic3.py


The SubplotWcs(or AxesWcs) allows you to plot in the sky
coordinates. For example, ``ax["fk5"]`` gives you an Axes whose data
coordinate is in fk5 coordinate. Most (if not all) of the valid mpl
plot commands will work. The unit for the sky coordinates are
degrees.::

  # (alpha, delta) in degree
   ax["fk4"].plot([x/24.*360 for x in [4, 5, 6]],
                  [-74, -70, -66], "ro-")

  # (l, b)  in degree
  ax["gal"].plot([(285), (276.)],
                 [(-30), (-36)])

.. plot:: figures/demo_basic4.py


Instead of string ("fk4", "fk5", "gal"), you can use other pywcs.WCS
instance. The returning axes has a data coordinate of the image
cooridnate of the given wcs.

Displaying images in other wcs coordinate system is a bit tricky. You
may simply use imshow, which will regrid the original image into the
target wcs. Vector drawing using pcolormesh is recommedned, but only
when you use agg backend (or pdf where pcolormesh will be rasterized
with agg backend). Otherwise (e.g., ps), it will be extremely slow to
draw. Contour will be drawn in the original wcs coordinate and then
will be transformed to the target coordinate.

Here is a more sophiscated example. The two images are plotted using
the mpl_toolkits.AxesGrid. Both axes are created using the wcs
information of the first image. Note that the gridhelper object is
explicitly created and handed to the axes, i.e., the gridhelper is
shared between two axes. This is  The second image, which has different wcs
infromation is drawn using pcolormesh.


.. plot:: figures/demo_skyview.py
   :include-source:

