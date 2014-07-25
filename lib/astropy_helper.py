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

    def pix2world(self,*args,**kwargs):
        return self.wcs_pix2sky(*args,**kwargs)

    def world2pix(self,*args,**kwargs):
        return self.wcs_sky2pix(*args,**kwargs)

    pywcs.WCS.wcs_pix2world = pix2world
    pywcs.WCS.wcs_world2pix = world2pix
