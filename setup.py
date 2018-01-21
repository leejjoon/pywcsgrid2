# from ez_setup import use_setuptools
# use_setuptools()

# from setuptools import setup
from setuptools import setup

def main():

    setup(name = "pywcsgrid2",
          version = "1.0-git",
          description = "pywcsgrid2 is a python module to be used with matplotlib for displaying astronomical fits images",
          author = "Jae-Joon Lee",
          author_email = "lee.j.joon@gmail.com",
          url="http://leejjoon.github.com/pywcsgrid2/",
          license = "MIT",
          platforms = ["Linux","Mac OS X"],
          packages = ['pywcsgrid2'],
          package_dir={'pywcsgrid2':'lib',
                       },
          use_2to3 = True,
          )


if __name__ == "__main__":
    main()
