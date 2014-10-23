"""
Aerosol opacity structure (vertical distribution) and scattering 
phase functions as measured by DISR on the Huygens probe.
"""

import numpy as np
import os

pi = 3.141592653589793

datapath = os.path.join(os.getenv('RTDATAPATH'),'aerosol/titan/')

try:
    Tomasko = {'lo':np.genfromtxt(datapath+'Tomasko2007_phase_0-80km.TAB',  names=True),
               'hi':np.genfromtxt(datapath+'Tomasko2007_phase_80-200km.TAB', names=True)}
except:
    None

def set_opacity(model, method='altitude'):
    """Set aerosol haze opacity for atmosphere data structure, model, using
    a verticual distribution that is consistent with the total Huygens/DISR
    opacity and either uniform in 'altitude' or 'pressure', as specified by 
    the method keyword."""
    
    if 'wavelength' not in model.keys():
        print('WARNING: Set wavelength scale first (e.g., with CH4 opacity.)')
        return None
    
    tau = {'low':6.270e2  * (1e3*model['wavelength'])**(-0.9706), 
           'mid':2.029E4  * (1e3*model['wavelength'])**(-1.409 ),     
           'high':1.012E7 * (1e3*model['wavelength'])**(-2.339 ),
           'method':method,}

    tau.update({'total':tau['low']+tau['mid']+tau['high'] })
    model.update({'haze':{'tau':tau}})
    
    i_low = np.where( model['layers']['z'] < 30)[0]
    i_mid = np.where((model['layers']['z'] >= 30) &
                     (model['layers']['z_max'] < 80) )[0]
    i_high = np.where(model['layers']['z_max'] >= 80)[0]

    tau_haze = np.empty((model['nlay'],model['nlam']))
        
    if method == 'pressure': # Uniform vertical opacity distribution with pressure
        tau_haze[i_low,:] = np.array([model['haze']['tau']['low']/len(i_low)]*len(i_low))
        tau_haze[i_mid,:] = np.array([model['haze']['tau']['mid']/len(i_mid)]*len(i_mid))
        tau_haze[i_high,:] = np.array([model['haze']['tau']['high']/len(i_high)]*len(i_high))
    
    if method == 'altitude': # Uniform vertical opacity distribution with altitude
        for i in i_low:
            tau_haze[i,:] = tau['low']*model['layers']['dz'][i]/30.0

        for i in i_mid:
            tau_haze[i,:] = tau['mid']*model['layers']['dz'][i]/50.0

        tau_haze[i_high,:] = np.array([tau['high']/len(i_high)]*len(i_high))
 
    model['layers']['tau'].update({'haze':tau_haze})
   
    return

def show__opacities_DISR(model, ax=None):
    """Show the three regions of aerosol opacity defined in
    Tomasko et al., 2008."""

    import matplotlib.pyplot as plt

    if not ax : fig, ax = plt.subplots()
    ax.plot(model['wavelength'], model['haze']['tau']['low'], 'k--')
    ax.plot(model['wavelength'], model['haze']['tau']['mid'], 'k:')
    ax.plot(model['wavelength'], model['haze']['tau']['high'], 'k-.')
    ax.plot(model['wavelength'], model['haze']['tau']['total'], 'k-', lw=2)
    ax.set_ylabel(r'aerosol optical dpeth') ; 
    ax.set_xlabel(r'wavelength ($\mu {\rm m}$)') ; 

def legendre_eval(x,coeff):
    """Evaluate a Legendre polynomial with coefficients,  
    coeff, while fixing the zero-th order coefficient at unity.

    Uses the recurrence relation of Legendre polynomials
        (n+1)*P_n+1(x) = (2n+1)*x*P_n(x) - n*P_n-1(x)
    evaluated with the Clenshaw recurrence formula, see Numerical Recipes
    by Press et al. (1992), Section 5.5
    
    Converted from polyleg.pro IDL procedure.
    """

    N  = len(coeff)-1 ; M = len(x)
    y  = np.zeros( (M,N+2) )
    jj = np.arange(N)+2
    
    coeff[0] = 1.0
    beta1 = -jj / (jj+1)
    for j in range(N,0,-1):
        alpha = (2*j + 1.)*x/(j + 1.)
        y[:,j-1] = (alpha*y[:,j] + beta1[j-1]*y[:,j+1] + coeff[j])
        
    return -0.5*y[:,1] + x*y[:,0] + coeff[0]

def phase_fit(deg, phase, order, efunc=None):
    """Fit Legendre polynomials to DISR phase function.
    
        phase   - scattering phase amplitude 
        deg     - scattering angle (degrees)
        order - maximum order of Legendre polynomial
        efunc   - phase functin error estimate
    """
    
    from scipy.optimize import leastsq
    from scipy.interpolate import interp1d

    residuals = lambda p, y, x, err: (y - legendre_eval(x,p))/err
    
    dom = np.cos(deg*pi/180)
    phase_interp = interp1d(deg, phase, bounds_error=False)

    degint = np.linspace(0,180,100)
    xint = np.cos(degint*pi/180.)
    yint = phase_interp(degint)

    if not efunc: efunc = lambda x: 0.10
    error = efunc(degint)*yint
    
    fit_coeffs = leastsq(residuals, np.ones(order+1), args=(yint, xint, error))
    
    return fit_coeffs[0]

def fit_DISR_phases():
    """Calculate Legendre fit coefficients and phase function 
    moments for Tomasko et al., scattering phase functions. 
    Return data structure with tabulated coefficients."""

    DISR = {'lambda':np.array([int(x) for x in Tomasko['lo'].dtype.names[1:]])/1e4}
 
    for alt in ['lo','hi']:
        phase_deg = Tomasko[alt]['deg']
        DISR.update({alt:{}})
        for nmom in [24,32]:
            coeff_table = np.empty((len(DISR['lambda']),nmom+1))
            pmom_table  = np.empty((len(DISR['lambda']),nmom+1))
            for i,wavelength in enumerate(DISR['lambda']):
                phase_amp =  Tomasko[alt]["{:.0f}".format(wavelength*1e4)]
                coeff_leg = phase_fit(phase_deg, phase_amp, nmom)
                pmom = coeff_leg/(2*np.arange(nmom+1)+1.)
                coeff_table[i,:] = coeff_leg
                pmom_table[i,:] = pmom
            DISR[alt].update({nmom:{'leg':coeff_table,
                                    'pmom':pmom_table,}})
    return DISR 


def show_phase_fit(wavelength, order=None, axs=None, efunc=None):
    """Plot DISR phase functions and Legendre fit at given wavelength. 
    Default order=64."""

    from matplotlib import gridspec
    
    if not axs: 
        fig = plt.figure(figsize=(5,6))
        subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1]) 
        ax = plt.subplot(gs[0])
        axr = plt.subplot(gs[1])
    else:
        ax = axs[0]
        axr = axs[1]
        
    order = 64
    if wavelength > 1000 : order = 48
    if wavelength > 1500 : order = 32        
    if wavelength > 3000 : order = 24
    if wavelength > 4000 : order = 16    
 
    if not efunc: efunc = lambda x: np.clip(x, 0.1, 0.1) 

    deg = Tomasko['lo']['deg'] 
    dom = np.cos(deg*3.141592653589793/180)
    phase= Tomasko['lo']["{:.0f}".format(wavelength*10.)]

    coeff = phase_fit(deg, phase, order, 
                      efunc=efunc)
    fit = legendre_eval(dom, coeff)
    ax.semilogy(deg, fit, 'k--', lw=2, label="fit <80 km, n={:d}".format(order))
    axr.plot(deg, (fit-phase)/phase, 'k--', lw=2,  drawstyle='steps-mid')
    axr.plot(deg, efunc(deg), 'k--', lw=1, color='gray')
    axr.plot(deg, -efunc(deg), 'k--', lw=1, color='gray')

    ax.scatter(deg, Tomasko['hi']["{:.0f}".format(wavelength*10.)], label='DISR >80 km', 
            s=80, marker='s', facecolor='none', edgecolor='black', linewidth=1)

    ax.scatter(deg, Tomasko['lo']["{:.0f}".format(wavelength*10.)], label='DISR <80 km',  
            s=80, marker='x', facecolor='none', edgecolor='black', linewidth=1)

    deg = Tomasko['hi']['deg'] 
    dom = np.cos(deg*3.141592653589793/180)
    phase= Tomasko['hi']["{:.0f}".format(wavelength*10.)]

    coeff = phase_fit(deg, phase, order, efunc=efunc)
    fit = legendre_eval(dom, coeff)
    ax.semilogy(deg, fit, 'k:', lw=2, label="fit >80 km, n={:d}".format(order))
    axr.plot(deg, (fit-phase)/phase, 'k:', drawstyle='steps-mid', lw=2, )

    ax.set_xlim(0,180)
    ax.text(10, 400, "{:4.1f}nm".format(wavelength), fontsize=12)

    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels(['','','-1','0','1','2','3'])
    ax.set_ylabel(r'phase function, log $P\,(\theta)$')
    ax.legend(fontsize=10, 
              scatterpoints=1,
              handlelength=3)
    ax.set_ylim(7e-2,1e3)
    axr.set_ylabel('residuals') ; axr.set_ylim(-0.25,0.25)
    axr.set_xlabel(r'scattering phase angle, $\theta$ (deg)')

def show_phase_function_fits(plotname=None):
    """Create figure with phase function fits at 
    at the 12 wavelengths between 0.7 -- 3.7um."""

    from matplotlib import gridspec
    fig = plt.figure(figsize=(18,15))
    gs0 = gridspec.GridSpec(3,4,wspace=0.3)

    for i in range(12):
        subgs = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[4,1], 
                                                 subplot_spec=gs0[i], hspace=0)

        axs = [plt.subplot(subgs[0]), plt.subplot(subgs[1]),]
        wav = int(Tomasko['lo'].dtype.names[i+5])/10
        show_phase_fit(wav, axs=axs,
                       efunc=lambda x: np.clip(x*2e-3, 0.1, 0.1))

    if plotname: savefig(plotname, bbox_inches='tight')

def phase_fit_coef_table_latex(DISR, nmom, alt, latex_file=None):
    """Create latex table of Legendre fit coefficients in 1-3um region."""

    tbl_txt = r"""
\begin{table}[ht]
\caption{Legendre coefficients fit to DISR phase functions}
\begin{center}
%\texttt{
\begin{tabular}{*{15}{r}} 
\hline\hline
"""

    wmn, wmx = 4,18    
    wavelengths = DISR['lambda'][wmn:wmx]
    
    tbl_txt += "\multicolumn{15}{c}{Wavelength ($\mu$m)} \\\\ \n"

    line = "{:^4s}".format('$k$')
    for wav in wavelengths:
        line += r"& {:^5.3f}".format(wav)
    line += "\\\\ \n \\hline \n"
    tbl_txt += line
    
    for k in range(nmom+1):
        line = "{:^4d}".format(k)
        for w in range(wmn,wmx):
            line += "&{:>7.4f} ".format(DISR[alt][nmom]['pmom'][w,k])
        line += '\\\\ \n'
        tbl_txt += line

    tbl_txt += """\hline\n\end{tabular}\n\n\end{center}\end{table}"""    

    if latex_file:
        with open(latex_file, "w") as f:
            f.write(tbl_txt)
        
    return tbl_txt


def set_aerosol_phase_moments(model, DISR, nmom):
    """Append aerosol scattering phase function moments,
    which are used by DISORT, to the atmospheric data
    structure. Coefficients are fit to DISR phase 
    functions and linearly interpolated between the 
    wavelengths that are tabulated in the literature."""

    from scipy.interpolate import interp1d

    wavelength = model['wavelength']
    
    pmom = {'lo':np.ones((len(wavelength),nmom+1)),
            'hi':np.ones((len(wavelength),nmom+1))}

    for alt in ['lo','hi']:
        for order in range(1,nmom+1) :
            f = interp1d(DISR['lambda'], 
                         DISR[alt][nmom]['pmom'][:,order],)
            pmom[alt][:,order] = f(wavelength)
       
    model['haze'].update({'pmom':pmom})


def show_aerosol_phase_moments(model, ax=None):
    """Plot DISR phase function moments together with
    linear interpolation."""
    fig, ax = subplots()
    alt = 'lo'
    nmom = len(model['haze']['pmom']['lo'][0,:])-1
    pmom = model['haze']['pmom']

    for order in range(1,8) :
        ax.plot(model['wavelength'], pmom['lo'][:,order], ls='--')
        ax.plot(DISR['lambda'], DISR[alt][nmom]['pmom'][:,order],
                'ko', alpha=0.5)

    ax.set_xlabel('wavelength (um')
    ax.set_ylabel('pmom')
    ax.set_xlim(min(model['wavelength']),
                max(model['wavelength']))    