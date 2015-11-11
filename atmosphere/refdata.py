"""
List of reference data files used by radiative transfer model to 
create a near-IR spectrum of Titan and methods to download data
from web and create reference data direcotry structure. 
"""

file_list = [
    'atmosphere_structure/titan/HASI_L4_ATMO_PROFILE_DESCEN.TAB',
    'aerosol/titan/Tomasko2007_phase_0-80km.TAB',
    'aerosol/titan/Tomasko2007_phase_80-200km.TAB',
    'gas_opacity/CIA/N2_N2.fits',
    'gas_opacity/CIA/H2_N2.fits',
    'gas_opacity/kc_CH4.SINFONI.v08.dnu_3.0.fits',
    'gas_opacity/kc_CH4.VIMS.v08.fits',
    ]


def get_from_url(fname, url='http://astro.berkeley.edu/~madamkov/refdata/'):
    """Download reference data from server."""

    import sys
    import os
    import urllib
    refdir = os.getenv('RTDATAPATH')
    ref_path = os.path.join(refdir,fname)

    if sys.version_info.major > 2:        
        try:
            print('Downloading...: {:s}'.format(url+fname))
            urllib.request.urlretrieve(url+fname, ref_path)
            print('          done.')
        except:
            print('Problem with {:s}'.format(url+fname))
    else:
        try:
            print('Downloading...: {:s}'.format(url+fname))
            urllib.urlretrieve(url+fname, ref_path)
            print('          done.')
        except:
            print('Problem with {:s}'.format(url+fname))
    return

def setup_directory():
    """Create reference directory structure, checking for exisiting
    data, and download from web if necessary"""

    import os
    refdir = os.getenv('RTDATAPATH') 

    for fname in file_list:
        ref_path = os.path.join(refdir,fname)
        ref_subdir = os.path.dirname(ref_path)
        if not os.path.exists(ref_subdir):
            os.makedirs(ref_subdir)
            print('Creating directory: {:s}'.format(ref_subdir))
        if not os.path.exists(ref_path): 
            get_from_url(fname)
        else:
            print('Reference data exists: {:s}'.format(ref_path))