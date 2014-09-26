"""
The parallel module doesn't work quite right yet. The namespaces 
associated with the engines in the notebook, for example, to know 
which engine 'id' is running in order to divide for loops into chunks,
is not available to the methods in the modules. However, this is not the 
case for methods defined in the notebook. Some other parallelization
scheme will need to be invoked
"""


def twostream_parallel(model, hg, mu0, rsurf):
    """Run two-stream calculation in parallel for
    atmosphere model and apend spectrum to data
    structure. Assumes that multiple engines are
    running in IPython notebook."""
    
    def calc_sub_spectrum(model, hg, mu0, rsurf):
        """Run two-stream calculation on a subset of the
        total spectral bandpass."""

        import numpy as np
        from atmosphere.rt import twostream    

        tau, ssa = twostream.get_opacity_ssa(model, verbose=False)
        intensities = np.zeros([model['nlam'], model['layers']['kc']['ng']])
        g_weight = model['layers']['kc']['w']

        chunk = model['nlam']/stride
        imn = int(id*chunk)
        imx = int((id+1)*chunk)
        
        if np.size(mu0) == 1:
            MU0 = np.empty(model['nlay']) 
            MU0.fill(mu0)
        else:
            MU0 = mu0

        for iwn in range(imn,imx):
            for g in range(model['layers']['kc']['ng']):
                tau_mono = tau[:,iwn,g]
                ssa_mono = ssa[:,iwn,g]
                g_Imu = twostream.mono(tau_mono, ssa_mono, hg, rsurf, MU0)
                intensities[iwn,g] = g_Imu[0]

        sub_spectrum = np.sum(intensities*g_weight, axis=1)  

        return sub_spectrum

    spectrum = np.zeros(model['nlam'])
    spectra = dv.apply_async(calc_sub_spectrum, model, hg, mu0, rsurf)

    for sub_spectrum in spectra.get():
        spectrum += sub_spectrum

    model.update({'spectrum':spectrum,
                'calc':{'method':'twostream_parallel',
                        'hg':hg,
                        'mu0':mu0,
                        'rsurf':rsurf,
                        }
                })