"""
Add gas opacities to model based on the composition and vertical structure.
"""

import numpy as np

def interpolate_kc(p, T, kc, verbose=False):
    """Linearly interpolate k-coefficients at a particular 
    pressure and temperature, using the input k-coefficent grid, kc.
    The standard structure of k-coefficient data array is:
    
        [wavelengths,pressures,temperatures,g-nodes]
    
    where the g-node are the Legendre-Gauss quadrature nodes
    or "g-ordinate". Returned array of coefficients corresponds to:
        
        [wavelengths,g-nodes]
    """
    
    pressures = np.array(kc['pressures'])
    temperatures = np.array(kc['temperatures'])
    
    ind_p = np.where(pressures < p)
    ind_T = np.where(temperatures < T)

    i = (np.max(ind_p) if np.size(ind_p) else np.array(0)).clip(0,len(pressures)-2) 
    j = (np.max(ind_T) if np.size(ind_T) else np.array(0)).clip(0,len(temperatures)-2)

    L11 = np.log(kc['kc'][:,i,j,:])
    L12 = np.log(kc['kc'][:,i+1,j,:])
    L21 = np.log(kc['kc'][:,i,j+1,:])
    L22 = np.log(kc['kc'][:,i+1,j+1,:])

    L1T = L11 + (L12-L11)*(T-temperatures[j])/(temperatures[j+1]-temperatures[j])
    L2T = L21 + (L22-L21)*(T-temperatures[j])/(temperatures[j+1]-temperatures[j])
    LPT = L1T + (L2T-L1T)*((np.log(p)-np.log(pressures[i]))/
                           (np.log(pressures[i+1])-np.log(pressures[i])))

    kc_interp = np.exp(LPT)

    if verbose:
        ik, ig = 100, 2
        print("Input p ={:8.2e} mbar,  T={:5.1f}K".format(p,T))
        print("'bracketting' array values:\n {:8.2e} -- {:8.2e} mbar\n {:8.1f} --{:8.1f}K".format(
                pressures[i], pressures[i+1], temperatures[j], temperatures[j+1]))
        print("A grid of kc-values at i_wavlength={:d}, ig={:d}:".format(ik, ig))
        print(kc['kc'][ik,i:i+2,j:j+2,ig])
        print("interpolated value: {:12.3e}".format(kc_interp[ik,ig]))
    return kc_interp

def append_kc_to_layers(model, kc, species):
    """Set k-coefficients for each layer by interpolating
    to appropriate temperature and pressure and update the 
    data structure for the amtosphere."""
    
    kc_shape = (model['nlay'], model['nlam'], kc['ng'])
    model['layers'].update({'kc':{species:np.ndarray(kc_shape),
                                'ng':kc['ng'],
                                'g':kc['g'],
                                'w':kc['w'],
                                }})    
    for i in range(model['nlay']):
        model['layers']['kc'][species][i,:,:] = interpolate_kc(model['layers']['p'][i],
                                                             model['layers']['T'][i],
                                                             kc)
    return


# There are additional k-coefficients for C2H2, C2H6, and CO. 
# These  are currently calculated on the VIMS grid as they are applicable in 
# roughly the 2.7--3um wavelength region that is inacccessible from the ground.
#
# Here we revise the gas opacity in the model to include multiple k-coefficient files. 
#
# It is a reasonable estimate to sum k-coefficients after interpolating each onto the 
# same pressure and temperature, however, a minimal amount of error-checking should confirm
# that the same wavelength grid and g-ordinates are being used.
#
# Overview of revisions to the code:
#
#    (1) Generalization of the set_methane() to other species.
#    (2) Error-checking for wavelength (and P,T) grid
#    (3) back-compatibility for set_methane() method.
#    (4) Some thought to CH3D abundance variability.
#

def set_methane(model, kc_file, CH3D_scale=None, verbose=False):
    """Set methane opacities in atmosphere structure, model, by 
    interpolatating k-coefficents from the specied kc_file,
    using the temperatures and pressures for each layer. 
    """

    if CH3D_scale:
        if len(kc_file) != 2:
            print('two k-coefficient files needed for set_methane_opacity()')
            return None
        kc = np.load(kc_file[0]).item()
        kc_CH3D = np.load(kc_file[1]).item()
        kc['kc'] = kc['kc']+CH3D_scale*kc_CH3D['kc']

        model.update({'wavelength':kc['wavelength']['mu'], 
                    'nlam':kc['wavelength']['nlam'], })                     
        append_kc_to_layers(model, kc, 'CH4')
        tau_CH4 = model['layers']['kc']['CH4'] * np.reshape(model['layers']['N_CH4'], 
                                                           (model['nlay'],1,1))

        if 'tau' not in model['layers']: model['layers'].update({'tau':{}})
        model['layers']['tau'].update({'CH4':tau_CH4})
        return

    if kc_file.endswith('.npy'): 
        kc = np.load(kc_file).item()
        model.update({'wavelength':kc['wavelength']['mu'], 
                    'nlam':kc['wavelength']['nlam'], })                     
        append_kc_to_layers(model, kc, 'CH4')
        tau_CH4 = model['layers']['kc']['CH4'] * np.reshape(model['layers']['N_CH4'], 
                                                           (model['nlay'],1,1))

        if 'tau' not in model['layers']: model['layers'].update({'tau':{}})
        model['layers']['tau'].update({'CH4':tau_CH4})
        return

    if kc_file.endswith('.fits'):   
        import pyfits
        hdu = pyfits.open(kc_file)

        kc = {'kc': hdu[0].data,
              'pressures':hdu[2].data['pressures'],
              'temperatures':hdu[3].data['temperatures'],
              'g': hdu[4].data['g'],
              'w': hdu[5].data['w'],
              'ng': hdu[0].header['NG'],
              }

        model.update({'wavelength':hdu[1].data['wavelength'], 
                      'nlam':len(hdu[1].data['wavelength']), 
                      })                        

        hdu.close()
        append_kc_to_layers(model, kc, 'CH4')
        tau_CH4 = model['layers']['kc']['CH4'] * np.reshape(model['layers']['N_CH4'], 
                                                           (model['nlay'],1,1))

        if 'tau' not in model['layers']: model['layers'].update({'tau':{}})
        model['layers']['tau'].update({'CH4':tau_CH4})

    return
    

def print_atmosphere_details(model):
    print('model dictionary data structure:')
    for item in model.keys():
        print("{0:7s} -  type: {2} - shape: {1}".format(
            item, shape(model[item]), type(model[item])))

    print("\natmosphere['layers'] dictionary data structure:")
    for item in model['layers'].keys():
        print("{0:7s} -  type: {2} - shape: {1}".format(
            item, shape(model['layers'][item]), type(model['layers'][item])))

def set_cia(model, scale=4.0, show_figure=False):
    """Append collision-induced-absorption opacity for
    N2-N2 and H2-N2 (in the near-IR) to the atmosphere 
    data structure, model."""

    import pyfits
    import os

    fits = pyfits.open(os.path.join(os.getenv('RTDATAPATH'), 
                                    'gas_opacity/CIA/N2_N2.fits'))
    k_N2N2 = fits[0].data
    fits.close()

    fits = pyfits.open(os.path.join(os.getenv('RTDATAPATH'), 
                                    'gas_opacity/CIA/H2_N2.fits'))
    k_H2N2 = fits[0].data
    fits.close()

    if 'wavelength' not in model.keys():
        print('WARNING: Set wavelength scale first (e.g., with CH4 opacity.)')
        return None
    
    tau_H2N2 = np.empty((model['nlay'],model['nlam']))
    tau_N2N2 = np.empty((model['nlay'],model['nlam']))
    
    layers = model['layers']
    N0 = 2.686e19 # Loschmidt number
    
    for i in range(model['nlay']):
        k_H2N2_interp = np.interp(model['wavelength'],
                                  (1e4/k_H2N2[::-1,0]), scale*k_H2N2[::-1,1],)

        k_N2N2_interp = np.interp(model['wavelength'],
                                  (1e4/k_N2N2[::-1,0]), scale*k_N2N2[::-1,1],)

        tau_H2N2[i,:] = k_H2N2_interp*layers['kmamg'][i] * \
                        layers['n'][i]/N0 * \
                        layers['m_N2'][i]*layers['m_H2']
        tau_N2N2[i,:] = k_N2N2_interp*layers['kmamg'][i] * \
                        layers['n'][i]/N0 * \
                        layers['m_N2'][i]*layers['m_N2'][i]

    layers['tau'].update({'H2_N2':tau_H2N2,
                          'N2_N2':tau_N2N2,})        
    if show_figure:
        fig, ax = subplots(figsize=(16,4))
        ax.plot(k_H2N2[:,0],  k_H2N2[:,1], 'k', drawstyle='steps-mid')
        ax.set_xlabel('wavenumber (cm$^{-1}$)')
        ax.set_ylabel(r'km$^{-1}$ amagat$^{-2}$')
        ax.set_xlim(4000,5000)

        fig, ax = subplots(figsize=(16,4))
        ax.plot(k_N2N2[:,0], k_N2N2[:,1], 'k', drawstyle='steps-mid')
        ax.set_xlabel('wavenumber (cm$^{-1}$)')
        ax.set_ylabel(r'km$^{-1}$ amagat$^{-2}$')
        ax.set_xlim(4000,5000) ;