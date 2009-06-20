from numpy.distutils.core import setup
from numpy.distutils.extension import Extension


import sys
from os.path import join

try:
    import numpy
except ImportError:
    print "numpy must be installed to build pywcsgrid."
    print "ABORTING."
    sys.exit(1)

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()



from numpy.distutils.misc_util import Configuration



import os.path



TPMDIR = './coords/tpm'

def get_tpm_sources(TPMDIR):
    import glob
    tpm_sources = glob.glob(os.path.join(TPMDIR, "*.c"))
    return tpm_sources


WCSLIBDIR="./wcslib-4.3.3/"
WCSLIBVER="4.3"

WCSLIBNAME = "wcs-%s"%WCSLIBVER


def main():

    setup(name = "pywcsgrid2",
          version = "0.1b1",
          description = "",
          author = "Jae-Joon Lee",
          maintainer_email = "lee.j.joon@gmail.com",
          license = "BSD",
          platforms = ["Linux","Mac OS X"],
          packages = ['pywcsgrid2'],
          package_dir={'pywcsgrid2':'lib',
                       },

          )


if __name__ == "__main__":
    main()
