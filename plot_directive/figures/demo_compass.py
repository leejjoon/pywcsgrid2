import pyfits
import matplotlib.pyplot as plt
import pywcsgrid2
import matplotlib.patheffects as patheffects


# use path_effects
my_path_effects = [patheffects.withStroke(linewidth=3,
                                          foreground="w")]
pywcsgrid2.Axes.default_path_effects = my_path_effects



if 1:

    f = pyfits.open("radio_21cm.fits")
    d, h = f[0].data, f[0].header

    plt.figure(1, [4.5,4.])

    ax = pywcsgrid2.axes([0.15, 0.15, 0.8, 0.8], header=h)

    ax.imshow(d, origin="low",
              cmap=plt.cm.gist_heat_r, interpolation="nearest")

    ax.set_xlim(20, 52)
    ax.set_ylim(20, 52)

    # Figure title
    t = ax.add_inner_title("Figure 1", loc=2, path_effects=False)

    # compass
    ax.add_compass(loc=1)

    # Beam size
    # (major, minor) = 3, 4 in pixel, angle=20
    ax.add_beam_size(3, 4, 20, loc=3)

    # Size
    cdelt = abs(h["cdelt1"])
    ax.add_size_bar(0.5/cdelt, # 30' in in pixel
                    r"$30^{\prime}$",
                    loc=8)

    plt.show()

