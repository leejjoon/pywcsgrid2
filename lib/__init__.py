

import re

_wcs_key_pattern = re.compile(r'^(NAXIS|CD|CDELT|CRPIX|CRVAL|CTYPE|CROTA|LONGPOLE|LATPOLE|PV|DISTORT|OBJECT|BUNIT|EPOCH|EQUINOX|LTV|LTM|DTV|DTM)')

def filterwcs(h):
    """
    select wcs related cards
    """
    import pyfits
    l = [card for card in h.ascardlist() if _wcs_key_pattern.match(card.key)]
    return pyfits.Header(l)



