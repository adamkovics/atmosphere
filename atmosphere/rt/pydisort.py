"""
Wrapper for CDISORT
"""

def mono(mono_input):
    """Monochromatic call to the CDISORT module."""

    import numpy as np
    import ctypes
    from atmosphere.rt import cdisort 

    ds = cdisort.disort_state() ; dout = cdisort.disort_output()

    TRUE, FALSE    = 1, 0 
    ds.accur       = 0.
    ds.flag.ibcnd  = FALSE
    ds.flag.usrtau = FALSE
    ds.flag.usrang = TRUE
    ds.flag.lamber = TRUE
    ds.flag.planck = FALSE
    ds.flag.onlyfl = FALSE
    ds.flag.quiet  = TRUE
    ds.flag.spher  = FALSE
    ds.flag.old_intensity_correction = TRUE

    ds.bc.umu0   = float(mono_input['umu0'])
    ds.bc.phi0   = float(mono_input['phi0'])

    ds.nstr   = int(mono_input['nstr'])
    ds.nphase = ds.nstr
    ds.nlyr   = int(mono_input['nlyr'])
    ds.nmom   = int(mono_input['nmom'])
    ds.numu   = 1
    ds.nphi   = 1

    cdisort.c_disort_state_alloc(ds)

    ds.bc.fbeam  = 1.0
    ds.bc.fisot  = 0.0
    ds.bc.albedo = float(mono_input['albedo'])

    ssalb  = ctypes.cast( ds.ssalb.__int__(), ctypes.POINTER( ctypes.c_double ) )
    umu    = ctypes.cast( ds.umu.__int__(),   ctypes.POINTER( ctypes.c_double ) )
    phi    = ctypes.cast( ds.phi.__int__(),   ctypes.POINTER( ctypes.c_double ) )
    dtauc  = ctypes.cast( ds.dtauc.__int__(), ctypes.POINTER( ctypes.c_double ) )
    pmom   = ctypes.cast( ds.pmom.__int__(),  ctypes.POINTER( ctypes.c_double ) )
        
    for i in range(ds.nlyr) : dtauc[i] = float(mono_input['dtau'][i])
    for i in range(ds.nlyr) : ssalb[i] = float(mono_input['ssalb'][i])

    #Scattering phase function moments the same in all layers
    phase_moments = mono_input['pmom']
    for i in range(ds.nlyr*(ds.nmom+1)):
        pmom[i] = float(phase_moments[np.mod(i, ds.nmom+1)])
        
    umu[0]   = float(mono_input['umue'])
    phi[0]   = float(mono_input['phie'])
        
    cdisort.c_disort_out_alloc(ds, dout)
    cdisort.c_disort(ds, dout)

    flup = dout.rad.flup
    uu   = ctypes.cast( dout.uu.__int__(), ctypes.POINTER( ctypes.c_double ) )[0:ds.numu*ds.ntau*ds.nphi]
    u0u  = ctypes.cast( dout.uu.__int__(), ctypes.POINTER( ctypes.c_double ) )[0:ds.numu*ds.ntau*ds.nphi]

    cdisort.c_disort_state_free(ds)
    cdisort.c_disort_out_free(ds, dout)

    return flup, u0u[0]


def get_tau_ssa(model, verbose=False):
    """Return cumulative layer optical depths and 
    single scattering abledos from atmosphere data structure."""

    import numpy as np

    def redimension_over_ng(array, ng):
        """Increase array dimensions by for summing opacities that 
        have uniform probability distribution in bandpass."""

        output = np.reshape(np.tile(np.reshape(array, array.shape+(1,)), ng),
                            array.shape+(ng,)
                            )
        return output  

    tau = np.zeros((model['nlay'], model['nlam']))
    ng = model['layers']['kc']['ng']
    
    kc_species = []
    for item in model['layers']['tau']:
        if len(model['layers']['tau'][item].shape) == 2:
            if verbose: print('Adding {:s} to 2d sum'.format(item))
            tau += model['layers']['tau'][item]
        else:
            if verbose: print('Adding {:s} to 3d list'.format(item))
            kc_species.append(item)
            
    tau_total = redimension_over_ng(tau, ng)
    tau_haze  = redimension_over_ng(model['layers']['tau']['haze'], ng)

    if verbose: print('New tau shape:', tau_total.shape)
    
    for item in kc_species:
        if verbose: print('Adding {:s} to total'.format(item))
        tau_total += model['layers']['tau'][item]

    ssa = tau_haze / tau_total
    
    return tau_total, ssa