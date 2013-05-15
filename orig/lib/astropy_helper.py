try:
    import astropy
except ImportError:
    has_astropy = False
else:
    has_astropy = True

if has_astropy:
    from astropy.io import fits as pyfits
    import astropy.wcs as pywcs
else:
    import pyfits
    import pywcs

