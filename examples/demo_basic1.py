
import pyfits
import matplotlib.pyplot as plt


f = pyfits.open("data/lmc.fits")
d, h = f[0].data, f[0].header

plt.figure(1, [6,3.5])

plt.subplot(121)

plt.imshow(d, origin="low", vmin=0, vmax=2000,
           cmap=plt.cm.gray_r)

plt.title("Original MPL")

import pywcsgrid2

pywcsgrid2.subplot(122, header=h)

plt.imshow(d, origin="low", vmin=0, vmax=2000,
           cmap=plt.cm.gray_r)

plt.title("pywcsgrid2")

plt.subplots_adjust(left=0.1, right=0.95, top=0.95, wspace=0.5)

plt.show()

