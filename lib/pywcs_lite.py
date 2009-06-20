import pywcs
import  pywcs._pywcs as _pywcs
import numpy

class Wcs(pywcs.WCS):

    def __init__(self, header=None, fobj=None, key=' ', minerr=0.0, relax=False, naxis=2):
        if naxis != 2:
            raise ValueError("Only 2 axes are supported")
        self.naxis = naxis

        if header is None:
            wcsprm = _pywcs._Wcsprm(header=None, key=key,
                                    relax=relax, naxis=naxis)
            # Set some reasonable defaults.
            wcsprm.crpix = numpy.zeros((self.naxis,), numpy.double)
            wcsprm.crval = numpy.zeros((self.naxis,), numpy.double)
            wcsprm.ctype = ['RA---TAN', 'DEC--TAN']
            wcsprm.cd = numpy.array([[1.0, 0.0], [0.0, 1.0]], numpy.double)
            cpdis = (None, None)
            sip = None
        else:
            try:
                header_string = "".join([str(x) for x in header.ascardlist()])
                _wcsprm = _pywcs._Wcsprm(header=header_string, key=key,
                                        relax=relax, naxis=naxis)
                wcsprm = _wcsprm.wcssub([_pywcs.WCSSUB_LONGITUDE,
                                         _pywcs.WCSSUB_LATITUDE]) # subtract celestail axes.
            except _pywcs.NoWcsKeywordsFoundError:
                wcsprm = _pywcs._Wcsprm(header=None, key=key,
                                        relax=relax, naxis=naxis)

            cpdis = self._read_distortion_kw(header, fobj, key=key,dist='CPDIS', err=minerr)
            sip = self._read_sip_kw(header, key=key)

        pywcs.WCSBase.__init__(self, sip, cpdis, wcsprm)
        self.footprint = self.calcFootprint(header)



if __name__ == "__main__":
    import pyfits
    f= pyfits.open("/export/bulk/blackstone1/lee/chandra/w63/G82.2+5.3.21cm.fits")
    wcs = Wcs(header=f[0].header)
