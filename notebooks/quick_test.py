import numpy as np
import os
from time import time

import sys
required_version = (3,3)
current_version = sys.version_info

os.environ['RTDATAPATH'] = '/Users/mate/g/rt/data/refdata/'

import atmosphere as atm
from atmosphere.rt import pydisort

atm.refdata.setup_directory()

def create_atmosphere_model(**kw):
    """
    Setup layered atmospheric model using a default physical/thermal structure
    determined by HASI, and composition from the GCMS, with aerosol 
    scattering propertier from DISR.
    
    The titan dictionary contains the values used to determine the opacity
    and scattering albedo of each atmospheric layer as a function of
    wavelength. The dispersion axis is constrained by the k-coefficients that
    are used, as specified by the methane kc_file.
    """
    
    ts = time()
    titan = atm.structure.set_HASI_structure(nlev=21,method='split_at_tropopause')
    atm.composition.set_abundances(titan, trace_gas={'m_H2':0.001})
    atm.gas_opacity.set_methane(titan, 
                                kc_file=os.path.join(os.environ['RTDATAPATH'],
                                                     'gas_opacity/kc_CH4.VIMS.v08.fits'),
                                )
    atm.gas_opacity.set_cia(titan)
    atm.aerosol.set_opacity(titan)
    DISR = atm.aerosol.fit_DISR_phases()
    atm.aerosol.set_aerosol_phase_moments(titan, DISR, nmom=32)   
    t_setup = time()-ts
    
    titan['haze'].update({'ssalb':0.96})
    titan.update({'radius':2575.,
                  'rsurf':0.10,
                  't_setup':t_setup,
                 })
    if 'rsurf' in kw: titan.update({'rsurf':kw['rsurf']})
    if 'verbose' in kw and kw['verbose']:
        pstr = 'PyDISORT Titan atmosphere structure and opacity setup: {:4.1f} sec'
        print(pstr.format(titan['t_setup']))    
    
    return titan

def setup_VIMS_calc(**kw):
    """Specify the viewing geometry and wavelegth range for the
    radiative transfer calculation."""

    titan = create_atmosphere_model(**kw)
    titan.update({'rt':{'spher':False,
                        'wav_range':(2.0,2.40),
                        'view':{'umu0':0.99,'umue': 0.90,'phi0': 10.0,'phie': 11.0},
                       }})
    for k in ['view','wav_range']: 
        if k in kw: titan['rt'].update({k:kw[k]})
            
    fi = lambda array, v: abs(array-v).argmin()
    wav_indices = lambda array, wavs: tuple([fi(array, mu) for mu in wavs])
    wav_mn, wav_mx = wav_indices(titan['wavelength'], titan['rt']['wav_range'])
    titan['rt'].update({'wav_mn':wav_mn,
                        'wav_mx':wav_mx,
                        'nlam':wav_mx-wav_mn+1,
                        })          
    return titan

if __name__ == '__main__':

    if current_version >= required_version:
        ts = time()
        print('PyDISORT Running quick test.')
        titan = setup_VIMS_calc(wav_range=(2.0,2.01), verbose=True ) 
        pydisort.calculate_spectrum(titan)
        print("PyDISORT Calculation for {:d} channels complete: {:4.1f} sec [success]".format(
            titan['rt']['nlam'],
            time()-ts),)
    else:
        print('PyDISORT Currently requires Python 3.3 or newer.')
