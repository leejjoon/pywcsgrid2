from distutils.core import setup


def main():

    setup(name = "pywcsgrid2",
          version = "0.1a2",
          description = "pywcsgrid2 is a python module to be used with matplotlib for displaying astronomical fits images",
          author = "Jae-Joon Lee",
          author_email = "lee.j.joon@gmail.com",
          url="http://leejjoon.github.com/pywcsgrid2/",
          license = "MIT",
          platforms = ["Linux","Mac OS X"],
          packages = ['pywcsgrid2','pywcsgrid2.axes_grid'],
          package_dir={'pywcsgrid2':'lib',
                       'pywcsgrid2.axes_grid':'axes_grid',
                       },

          )


if __name__ == "__main__":
    main()
