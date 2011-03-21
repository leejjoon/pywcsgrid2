import matplotlib.pyplot as plt

from pywcsgrid2.allsky_axes import make_allsky_axes

if 1:

    # proj must be one of "CYP", "CEA", "CAR", "MER", "SFL", "PAR", "MOL"

    fig = plt.figure()

    rect = 211
    proj = "PAR"
    coord, lon_center = "fk5", 180
    ax1 = make_allsky_axes(fig, rect, coord=coord, proj=proj, lon_center=lon_center)
    ax1.grid()

    rect = 212
    proj = "CAR"
    coord, lon_center = "gal", 0
    ax2 = make_allsky_axes(fig, rect, coord=coord, proj=proj, lon_center=lon_center)
    ax2.grid()
    
    plt.show()
