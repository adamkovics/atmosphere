import numpy as np
from scipy.linalg import solve_banded

def mono(tau, ssa, hg, rsurf, MU0, F0PI=1.0, BTOP=0.0):
    """Monochromatic two-stream radiative transfer solver with:
           tau - optical depths per layer
           ssa - single scattering albedo
           hg  - scattering asymetry parameter
           rsurf - surface reflectivity
           MU0 - cos(theta) where theta is both the 
                 incidence an emission angle.
           F0PI - top of atmosphere incident flux * pi
           BTOP - top of the atmosphere diffuse flux.
    """

    import numpy as np
    nlay = len(tau)

    # Cumulative optical depth
    taub = np.cumsum(tau)
    taut = np.append(0,taub[:-1])

    # Surface reflectance and lower boundary condition
    bsurf = rsurf * MU0[-1] * F0PI * np.exp(-np.sum(tau)/MU0[-1]) 

    twopi = np.pi+np.pi
    g1 = 0.86602540378 * (2.-ssa*(1+hg))  
    g2 = (1.7320508075688772*ssa/2.) * (1-hg)       
    g3 = (1.-1.7320508075688772*hg*MU0)/2        
    g4 = 1. - g3   
    
    lam = np.sqrt(g1*g1 - g2*g2)
    gamma = (g1-lam)/g2
    alpha = np.sqrt( (1.-ssa) / (1.-ssa*hg) )   

    Am = F0PI * ssa *(g4 * (g1 + 1./MU0) + g2*g3 )/ (lam*lam - 1./(MU0*MU0))
    Ap = F0PI * ssa *(g3 * (g1 - 1./MU0) + g2*g4 )/ (lam*lam - 1./(MU0*MU0))

    # Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer.
    Cpm1 = Ap * np.exp(-taut/MU0) 
    Cmm1 = Am * np.exp(-taut/MU0) 
    # Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    Cp = Ap * np.exp(-taub/MU0) 
    Cm = Am * np.exp(-taub/MU0) 
        
    #  Solve for the coefficients of system of equations using boundary conditions
    # Exponential terms: 
    exptrm = lam*tau
    exptrm[exptrm>35] = 35 # clipped so that exponential doesn't explode
    Ep = np.exp(exptrm)
    Em = 1./Ep
    
    E1 = Ep + gamma*Em
    E2 = Ep - gamma*Em
    E3 = gamma*Ep + Em
    E4 = gamma*Ep - Em
    
    L = nlay+nlay
    Af = np.empty(L) 
    Bf = np.empty(L) 
    Cf = np.empty(L) 
    Df = np.empty(L)
 
    # First Term
    Af[0] = 0.0            
    Bf[0] = gamma[0] + 1.
    Cf[0] = gamma[0] - 1.
    Df[0] = BTOP - Cmm1[0]
    
    AA = (E1[:-1]+E3[:-1])*(gamma[1:]-1)
    BB = (E2[:-1]+E4[:-1])*(gamma[1:]-1)
    CC = 2.*(1.-gamma[1:]*gamma[1:])
    DD = (gamma[1:]-1) * (Cpm1[1:] - Cp[:-1]) + (1-gamma[1:]) * (Cm[:-1]-Cmm1[1:])
    Af[1:-1:2]=AA
    Bf[1:-1:2]=BB
    Cf[1:-1:2]=CC
    Df[1:-1:2]=DD

    AA = 2.*(1.-gamma[:-1]*gamma[:-1])
    BB = (E1[:-1]-E3[:-1])*(gamma[1:]+1.)
    CC = (E1[:-1]+E3[:-1])*(gamma[1:]-1.)
    DD = E3[:-1]*(Cpm1[1:] - Cp[:-1]) + E1[:-1]*(Cm[:-1] - Cmm1[1:])
    Af[2::2]=AA
    Bf[2::2]=BB
    Cf[2::2]=CC
    Df[2::2]=DD
    # Last term:
    Af[-1] = E1[-1] - rsurf*E3[-1]
    Bf[-1] = E2[-1] - rsurf*E4[-1]
    Cf[-1] = 0.0
    Df[-1] = bsurf - Cp[-1] + rsurf*Cm[-1]
    
    # Match to IDL 
    Af2=np.append(Af[1:],0)
    Cf2=np.append(0,Cf[0:-1])

    # Solve tridiagonal matrix
    ab=np.matrix([
        Cf2,
        Bf,
        Af2,
    ])
    
    k = solve_banded((1,1),ab,Df)
     
    # Unmix coefficients
    even = np.arange(0,2*nlay,2)
    odd  = even+1    
    k1 = k[even] + k[odd] 
    k2 = k[even] - k[odd]        
    
    #  Iteratively determine upward intensity (along mu) in each layer
    Fsurf = k1[-1]*Ep[-1] + gamma[-1]*k2[-1]*Em[-1] + Cp[-1] # flux at surface
    Isurf = Fsurf/np.pi           # Lambert intensity distribution
      
    Am = Am * np.exp(-taut / MU0)
    Ap = Ap * np.exp(-taut / MU0)

    F1 = 1 + 1.5 * hg * MU0
    F2 = 1 - 1.5 * hg * MU0

    G = (ssa * k1 * (F1 + gamma*F2)) / (twopi)          
    H = (ssa * k2 * (gamma*F1 + F2)) / (twopi)
    A = (ssa * (F1*Ap + F2*Am)) / (twopi)

    # Integrate upward intensity starting from surface 
    Imu = np.empty(nlay+1)        # upward intensity at each level (not layer)
    Imu[-1] = Isurf               # upward intensity from surface
    
    addterms = (
        (ssa*F0PI/(8*np.pi))           
        *  (1-hg)*np.exp(-taut/MU0) * (1 - np.exp(-2.*tau/MU0))      
        +  A * (1 - np.exp(-2*tau/MU0) ) / 2                            
        +  G * (np.exp(exptrm-tau/MU0) - 1) / (lam*MU0 - 1) 
        +  H * (1 - np.exp(-exptrm-tau/MU0))/ (lam*MU0 + 1)
        )
    multterms = np.exp(-tau/MU0)
    for j in np.arange(nlay-1,-1,-1):
        Imu[j] = Imu[j+1]*multterms[j] + addterms[j]
                          
    return Imu

def calc_spectrum(model, hg, mu0, rsurf):
    """Run two-stream calculation on a subset of the
    total spectral bandpass."""

    import numpy as np

    tau, ssa = get_opacity_ssa(model, verbose=False)
    intensities = np.zeros([model['nlam'], model['layers']['kc']['ng']])
    g_weight = model['layers']['kc']['w']

    imn = 0
    imx = int(model['nlam'])
    
    if np.size(mu0) == 1:
        MU0 = np.empty(model['nlay']) 
        MU0.fill(mu0)
    else:
        MU0 = mu0

    for iwn in range(imn,imx):
        for g in range(model['layers']['kc']['ng']):
            tau_mono = tau[:,iwn,g]
            ssa_mono = ssa[:,iwn,g]
            g_Imu = mono(tau_mono, ssa_mono, hg, rsurf, MU0)
            intensities[iwn,g] = g_Imu[0]

    spectrum = np.sum(intensities*g_weight, axis=1)  

    model.update({'spectrum':spectrum,
                  'calc':{'method':'twostream',
                        'hg':hg,
                        'mu0':mu0,
                        'rsurf':rsurf,
                        },
                 })

    return model


def get_opacity_ssa(model, verbose=False):
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


def pseudo_pp(mu0, model, R_surface):
    """Return incidence angles for pseudo-plane-parallel
    correction for spherical geometry of an extended atmosphere.
    Analogous to curved path through parallel layers. 
    Solid surface radius, R_surface (in km) and the altitude
    above the surface from the atmosphere data structure, atm."""
        
    mu0_sphere = np.sqrt(1-(R_surface/(R_surface+model['layers']['z']))**2*(1.0-mu0**2))
    
    return mu0_sphere    
