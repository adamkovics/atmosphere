"""

Methods for running CDISORT radiative tranfer via Python/SWIG wrappers
for atmospheric data structures generated with the atmosphere package.

"""

__version__ = '0.0.4'

from time import time
import numpy as np

def mono(mono_input):
    """Monochromatic call to the CDISORT module."""
    import numpy as np
    import ctypes
    from atmosphere.rt import cdisort 

    ds = cdisort.disort_state() 
    dout = cdisort.disort_output()

    TRUE, FALSE    = 1, 0 
    ds.accur       = 0.
    ds.flag.ibcnd  = FALSE
    ds.flag.usrtau = FALSE
    ds.flag.usrang = TRUE
    ds.flag.lamber = TRUE
    ds.flag.planck = FALSE
    ds.flag.onlyfl = FALSE
    ds.flag.quiet  = TRUE
    ds.flag.spher  = mono_input['spher']
    ds.flag.old_intensity_correction = TRUE
    ds.radius  = float(mono_input['radius'])
    ds.bc.umu0 = float(mono_input['umu0'])
    ds.bc.phi0 = float(mono_input['phi0'])
    ds.nstr   = int(mono_input['nstr'])
    ds.nphase = ds.nstr
    ds.nlyr   = int(mono_input['nlyr'])
    ds.nmom   = int(mono_input['nmom'])
    ds.numu   = 1
    ds.nphi   = 1

    cdisort.c_disort_state_alloc(ds)

    ds.bc.fbeam  = np.pi
    ds.bc.fisot  = 0.0
    ds.bc.albedo = float(mono_input['albedo'])

    ssalb  = ctypes.cast( ds.ssalb.__int__(), ctypes.POINTER( ctypes.c_double ) )
    umu    = ctypes.cast( ds.umu.__int__(),   ctypes.POINTER( ctypes.c_double ) )
    phi    = ctypes.cast( ds.phi.__int__(),   ctypes.POINTER( ctypes.c_double ) )
    dtauc  = ctypes.cast( ds.dtauc.__int__(), ctypes.POINTER( ctypes.c_double ) )
    pmom   = ctypes.cast( ds.pmom.__int__(),  ctypes.POINTER( ctypes.c_double ) )

    if mono_input['spher']:
        zd = ctypes.cast( ds.zd.__int__(),    ctypes.POINTER( ctypes.c_double ) )    
        for i in range(ds.nlyr) : zd[i] = float(mono_input['z'][i])
        zd[ds.nlyr] = 0.

    for i in range(ds.nlyr) : dtauc[i] = float(mono_input['dtau'][i])
    for i in range(ds.nlyr) : ssalb[i] = float(mono_input['ssalb'][i])
    for i in range(ds.nlyr*(ds.nmom+1)): pmom[i] = float(mono_input['pmom'][i])

    umu[0]   = float(mono_input['umue'])
    phi[0]   = float(mono_input['phie'])

    cdisort.c_disort_out_alloc(ds, dout)
    cdisort.c_disort(ds, dout)

    uu     = ctypes.cast(dout.uu.__int__(), 
                         ctypes.POINTER(ctypes.c_double))[0:ds.numu*ds.ntau*ds.nphi]
    u0u    = ctypes.cast(dout.u0u.__int__(), ctypes.POINTER(
                         ctypes.c_double))[0:ds.numu*ds.ntau*ds.nphi]
    mono_output = { 
        'uu':uu,         # Intensity, UU(iu,lu,j) (if ds.flag.onlyfl = FALSE; zero otherwise)
        'u0u':u0u,       # Az.av.Int. U0U(iu,lu,j) (if ds.flag.onlyfl = FALSE; zero otherwise)
        }

    cdisort.c_disort_state_free(ds)
    cdisort.c_disort_out_free(ds, dout)

    return mono_output

def get_tau_ssa(model, w_haze=0.96, verbose=False):
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

    ssa = w_haze * (tau_haze / tau_total)
    
    return tau_total, ssa

def get_disort_inputs(model) :
    """
    Generate the data structures for the monochromatic call
    to CDISORT. Assemble each channel into a list with 
    corresponding wavelength channel indices and g-ordinates,
    to be passed a multiprocessing pool.
    """
        
    from numpy import zeros, sum, abs, ndarray, floor, mod
    from atmosphere.rt.pydisort import get_tau_ssa, mono
    from copy import deepcopy

    mono_input = {'spher':model['rt']['spher'],
                  'radius':model['radius'],
                  'z':model['layers']['z_max'],
                  'nlyr':model['nlay'],
                  'umu0':model['rt']['view']['umu0'],
                  'phi0':model['rt']['view']['phi0'],
                  'umue':model['rt']['view']['umue'],
                  'phie':model['rt']['view']['phie'],
                  }
    i_mn = int(model['rt']['wav_mn'])
    i_mx = int(model['rt']['wav_mx'])

    tau, ssalb = get_tau_ssa(model, model['haze']['ssalb'])

    disort_mono_inputs = []
    for iwn in range(i_mn,i_mx):
        nmom = model['haze']['pmom']['hi'].shape[1]-1
        pmom = ndarray(model['nlay']*(nmom+1))
        for imom in range(model['nlay']*(nmom+1)):
            iz = int(floor(imom/(nmom+1)))
            if model['layers']['z'][iz] > 80 :
                phase_moments = model['haze']['pmom']['hi'][iwn,:]
            else:
                phase_moments = model['haze']['pmom']['lo'][iwn,:]
            pmom[imom] = float(phase_moments[mod(imom, nmom+1)])            
        mono_input.update({'nstr':nmom,
                           'nmom':nmom,
                           'albedo':model['rsurf'],
                           'pmom':pmom,
                           })
        for g in range(model['layers']['kc']['ng']) :
            mono_input.update({'dtau':tau[:,iwn,g],'ssalb':ssalb[:,iwn,g],})
            disort_mono_inputs.append((deepcopy(mono_input), iwn, g))

    return disort_mono_inputs

def calculate_spectrum(model) :
    """
    Loop over wavelength channels in model structure and
    call CDISORT via a single process.
    """
        
    from numpy import zeros, sum, abs, ndarray, floor, mod
    from atmosphere.rt.pydisort import get_tau_ssa, mono

    ts = time()
    mono_input = {'spher':model['rt']['spher'],
                  'radius':model['radius'],
                  'z':model['layers']['z_max'],
                  'nlyr':model['nlay'],
                  'umu0':model['rt']['view']['umu0'],
                  'phi0':model['rt']['view']['phi0'],
                  'umue':model['rt']['view']['umue'],
                  'phie':model['rt']['view']['phie'],
                  }
    
    intensity = zeros([model['nlam'], model['layers']['kc']['ng']]) 
    
    i_mn = int(model['rt']['wav_mn'])
    i_mx = int(model['rt']['wav_mx'])
    model['rt'].update({'parallel':'NOT implemented'})
    tau, ssalb = get_tau_ssa(model, model['haze']['ssalb'])
    for iwn in range(i_mn,i_mx):
        nmom = model['haze']['pmom']['hi'].shape[1]-1
        pmom = ndarray(model['nlay']*(nmom+1))

        for imom in range(model['nlay']*(nmom+1)):
            iz = int(floor(imom/(nmom+1)))
            if model['layers']['z'][iz] > 80 :
                phase_moments = model['haze']['pmom']['hi'][iwn,:]
            else:
                phase_moments = model['haze']['pmom']['lo'][iwn,:]
            pmom[imom] = float(phase_moments[mod(imom, nmom+1)])            
            
        mono_input.update({'nstr':nmom,
                           'nmom':nmom,
                           'albedo':model['rsurf'],
                           'pmom':pmom,
                           }
                          )
        for g in range(model['layers']['kc']['ng']) :
            mono_input.update({'dtau':tau[:,iwn,g],
                               'ssalb':ssalb[:,iwn,g],
                               })
            mono_output = mono(mono_input)
            intensity[iwn,g] = mono_output['uu'][0]

    spectrum = sum(intensity*model['layers']['kc']['w'], axis=1)
    model.update({'spectrum':spectrum})
    tc = time()-ts
    pstr = 'PyDISORT radiative transfer for {:d} channels: {:.1f} sec'       
    model['rt'].update({'method':'PyDISORT single process',
                        't_calc':tc,
                        'info':pstr.format(model['rt']['nlam'], 
                                           tc),
                        })
    return

def mp_spectrum_calculation(model, n_processes=4):
    """
    Calculate spectrum using the multiprocessing module
    with n_processes.
    """
    from multiprocessing import Pool
    from atmosphere.rt.pydisort import get_disort_inputs, mono

    ts = time()
    disort_inputs = get_disort_inputs(model)
    intensity = np.zeros([model['nlam'], model['layers']['kc']['ng']]) 

    with Pool(processes=n_processes) as pool:

        disort_outputs = [(pool.apply_async(mono, args=(channel,)),i,g) for 
                          channel,i,g in disort_inputs]

        for channel, i, g in disort_outputs:
            refl = channel.get(timeout=1) 
            intensity[i,g] = refl['uu'][0]

        spectrum = np.sum(intensity*model['layers']['kc']['w'], axis=1)
        model.update({'spectrum':spectrum})

    tc = time()-ts
    pstr = 'PyDISORT parallel RT ({:d} processes) for {:d} channels: {:.1f} sec'       
    model['rt'].update({'method':'PyDISORT parallel',
                        't_calc':tc,
                        'info':pstr.format(n_processes, 
                                           model['rt']['nlam'], 
                                           tc),
                        })    
    print("PyDISORT Calculation complete: {:7.1f}".format(time()-ts))

def ipcluster_spectrum_calculation(model):
    """
    Use the IPython parallel processing module to 
    run calculation on multiple engines.
    """
    from numpy import zeros
    from ipyparallel import Client, interactive

    def async_spectrum(model) :
        """
        This method is just like calculate_spectrum() above, but 
        the loop over wavelength has a range of wav_indices, which is
        distributed among engines using dv.scatter().
        """
        from numpy import zeros, sum, abs, ndarray, floor, mod
        from atmosphere.rt.pydisort import get_tau_ssa, mono

        mono_input = {'spher':model['rt']['spher'],
                      'radius':model['radius'],
                      'z':model['layers']['z_max'],
                      'nlyr':model['nlay'],
                      'umu0':model['rt']['view']['umu0'],
                      'phi0':model['rt']['view']['phi0'],
                      'umue':model['rt']['view']['umue'],
                      'phie':model['rt']['view']['phie'],
                      }

        intensity = zeros([model['nlam'], model['layers']['kc']['ng']]) 
        tau, ssalb = get_tau_ssa(model, model['haze']['ssalb'])
     
        for iwn in wav_indices:
            nmom = model['haze']['pmom']['hi'].shape[1]-1
            pmom = ndarray(model['nlay']*(nmom+1))

            for imom in range(model['nlay']*(nmom+1)):
                iz = int(floor(imom/(nmom+1)))
                if model['layers']['z'][iz] > 80 :
                    phase_moments = model['haze']['pmom']['hi'][iwn,:]
                else:
                    phase_moments = model['haze']['pmom']['lo'][iwn,:]
                pmom[imom] = float(phase_moments[mod(imom, nmom+1)])            
                
            mono_input.update({'nstr':nmom,
                               'nmom':nmom,
                               'albedo':model['rsurf'],
                               'pmom':pmom,
                               }
                              )
            for g in range(model['layers']['kc']['ng']) :
                mono_input.update({'dtau':tau[:,iwn,g],
                                   'ssalb':ssalb[:,iwn,g],
                                   })
                mono_output = mono(mono_input)
                intensity[iwn,g] = mono_output['uu'][0]

        return sum(intensity*model['layers']['kc']['w'], axis=1)

    ts = time()
    engines = Client()
    dv = engines[:]
    dv.scatter('wav_indices',range(model['rt']['wav_mn'],model['rt']['wav_mx']))

    spectra = dv.apply_async(interactive(async_spectrum), model)
    spectrum = zeros(model['nlam'])
    for sub_spectrum in spectra.get(): 
        spectrum += sub_spectrum

    model['rt'].update({'method':'parallel',
                        'parallel':'implemented',})
    model.update({'spectrum':spectrum})
    tc = time()-ts
    pstr = 'PyDISORT parallel RT ({:d} processes) for {:d} channels: {:.1f} sec'       
    model['rt'].update({'method':'PyDISORT ipyparallel',
                        't_calc':tc,
                        'info':pstr.format(len(dv), 
                                           model['rt']['nlam'], 
                                           tc),
                        })
    return 