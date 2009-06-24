========
Overview
========

The pywcsgrid2 is a collection of helper classes to display
astronomical fits images in matplotlib. Its main functionality is to
draw ticks, ticklabels, and grids in an appropriate sky
coordinates. I

.. contents::
   :depth: 1
   :local:


You require pyfits and kapteyn package installed, in addition to matplotlib.
You read in the fits file using pyfits and create a subplot (or axes) to
display the image with its header (wcs) information. ::

    import pyfits
    import pywcs
    import matplotlib.pyplot as plt
    from pywcsgrid2.axes_wcs import AxesWcs

    f = pyfits.open("data/lmc.fits")
    h = f[0].header

    fig = plt.figure(1, [5,5])
    ax = AxesWcs(fig, [0.2, 0.15, 0.7, 0.8], header=h)
    fig.add_axes(ax)

The data coordinate of the AxesWcs axes is the image (pixel) coordinate
((0-based!) of the fits file. But The ticklabels are displayed in the
sky coordinates. For a simple image display, you may do ::

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

Custom tick location is not currently (but wiil be) supported.

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


Again, the data coordinate of SubplotWcs(or AxesWcs) is a pixel
coordinate (0-based) of the fits header (or any equivalent wcs information).
For plot something in sky coordinate, you may convert your data into pixel coodinates, or you may use parasites axes which does that conversion for you. For
example, ``ax["fk5"]`` gives you an Axes whose data coordinate is in
fk5 coordinate. Most (if not all) of the valid mpl plot commands will
work. The unit for the sky coordinates are degrees.::

  # (alpha, delta) in degree
   ax["fk4"].plot([x/24.*360 for x in [4, 5, 6]],
                  [-74, -70, -66], "ro-")

  # (l, b)  in degree
  ax["gal"].plot([(285), (276.)],
                 [(-30), (-36)])

.. plot:: figures/demo_basic4.py


Instead of string ("fk4", "fk5", "gal"), you can use other pyfits
header instance. The returning axes has a data coordinate of the pixel
(image) cooridnate of the given header. Most of plot commands (other
than image-related routine) will work as expected.

However, displaying images in other wcs coordinate system needs some
consideration. You may simply use imshow:
   f2 = pyfits.open("another.fits")
   h2, d2 = f2[0].header, f2[0].data
   ax[h2].imshow(d2)

This will regrid the original image into the target wcs (regridin is
necessary since matplotlib's imshow only supports rectangular
image). If you don't want your data to be regridded, a vector drawing
command pcolormesh is recommedned. But pcolormesh is optimized for agg
backend and become extremely slow with increasing image size in other
backends. Therefore, it is highly recommended (and this is a default
behavior) that pcolormesh command is rasterized (rasterization is
fully supported in pdf and svg backend, and partially available in ps
backend). Contouring command will work fine. Contours will be drawn in
the original wcs coordinate and then will be transformed to the target
coordinate.

The code below is a more sophiscated example. The two fits images
are plotted using the mpl_toolkits.AxesGrid. Both axes are created
using the wcs information of the first image. Note that the gridhelper
object is explicitly created and handed to the axes, i.e., the
gridhelper is shared between two axes (this is to share grid
parameters). The second image, which has different wcs
infromation is drawn using pcolormesh.


.. plot:: figures/demo_skyview.py
   :include-source:

It is possible to create a floating axis in any sky coordinate. This
can be useful for drawing a Galactic object, where you draw a image in
RA-Dec, but want to indicate the Galactic location of the object. A
floating axis is created using the new_floating_axis method. The first
argument indicate which coordinate, and the second argument is the
value. For example, if you want to have a floating axis of b=0,
i.e. the second coordinate (index starts at 0) is 0 in the Galactic
coordinate:: 
  ax1.axis["b=0"] = ax["gal"].new_floating_axis(1, 0.)

Here is an complete example,

.. plot:: figures/demo_floating_axis.py
   :include-source:

