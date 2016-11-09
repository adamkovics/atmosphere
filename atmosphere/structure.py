import numpy as np
import os
datapath  = os.getenv('RTDATAPATH')
file_HASI = os.path.join(datapath,
                         'atmosphere_structure/titan/',
                         'HASI_L4_ATMO_PROFILE_DESCEN.TAB')

dtype_HASI = [('time', int), # time [milliseconds] 
              ('z', float),  # altitude [m]
              ('p', float),  # pressure [Pa]
              ('T', float),  # temperature [K]
              ('n', float)]  # densitgy [kg/m^3]
try:
    HASI = np.genfromtxt(file_HASI, dtype=dtype_HASI) 
except:
    url = 'http://w.astro.berkeley.edu/~madamkov/refdata/atmosphere_structure/titan/HASI_L4_ATMO_PROFILE_DESCEN.TAB'
    HASI = np.genfromtxt(url, dtype=dtype_HASI) 

HASI['z'] *= 1e-3 # convert [m] to [km]
HASI['p'] *= 0.01 # convert [Pa] to [mbar]
HASI['n'] *= 1e3 * 1e-6 * 6.0221e23 / 27.8  # convert [kg/m^3] to [cm-3]
        
amg = 2.6867805e19 # [cm-3]

def HASI_conversion(alt, outp, warnings=False):
    """Wrapper to create methods that interpolate the HASI atmospheric 
    structure specified by 'outp' (i.e., 'p','z','T','n') at an altitude 
    specified by 'p' or 'z'.Extropolation outside the HASI bounds uses 
    the boundary values.
    """ 
    from scipy.interpolate import interp1d
    if alt  not in ['p','z']: raise ValueError(alt)
    if outp not in ['p','z','T','n' ]: raise ValueError(outp)

    def conv_func(z):
        func = interp1d(HASI[alt], HASI[outp])
        try:
            return func(z)
        except ValueError:
            if warnings: print('WARNING: z outside of HASI range, extrapolating.')
            func = interp1d(HASI[alt], HASI[outp], bounds_error=False)
            out = func(z)
            if out.size == 1:
                if alt == 'p':
                    if z > HASI[alt][-1]: out = HASI[outp][-1]
                    if z < HASI[alt][0]:  out = HASI[outp][0]
                else:
                    if z < HASI[alt][-1]: out = HASI[outp][-1]
                    if z > HASI[alt][0]:  out = HASI[outp][0]
            else:
                if alt == 'p':
                    out[np.where(np.isnan(out) & (np.array(z) > HASI[alt][-1]))] = HASI[outp][-1]
                    out[np.where(np.isnan(out) & (np.array(z) < HASI[alt][0]))] = HASI[outp][0]
                else:
                    out[np.where(np.isnan(out) & (np.array(z) < HASI[alt][-1]))] = HASI[outp][-1]
                    out[np.where(np.isnan(out) & (np.array(z) > HASI[alt][0]))] = HASI[outp][0]
            return(out)
 
    alt_str = {'p':'pressure, z [mbar]',
               'z':'altitude, z [km]'}
    outp_str = {'p':'pressure [mbar]',
                'z':'altitude [km]',
                'T':'temperature [K]',
                'n':'number density[cm-3]',}

    conv_func.__doc__ = """The {:s} on Titan at {:s}, via 
    linear interpolation of HASI atmospheric structure.
    Output array dimensions are the same as input, z. 
    Values for z outside of HASI range return {:s} at
    boundary values.""".format(outp_str[outp], alt_str[alt], outp)
    
    return conv_func

temperature_at_pressure = HASI_conversion('p', 'T')
temperature_at_altitude = HASI_conversion('z', 'T')
density_at_pressure = HASI_conversion('p', 'n')
density_at_altitude = HASI_conversion('z', 'n')
pressure_at_altitude = HASI_conversion('z', 'p')
altitude_at_pressure = HASI_conversion('p', 'z')


def fi(array, v): 
    """Find array index (or indices) nearest to value(s) of v.
       Returns np.int64 or list."""
    if type(v) == float or type(v) == int: 
        return (np.abs(array-v)).argmin()
    else:
        i = []
        for vi in v: i.append((np.abs(array-vi)).argmin())
        return i

def set_HASI_structure(nlev=25, method='split_at_tropopause'):
    """Define layers in atmosphere based on HASI atmospheric structure,
    using one of the following methods:    
    
    linear_pressure | log_pressure | split_at_tropopause | equal_columns  

    Values (units) returned for each layer:
        p - pressure (mbar)
        T - temperature (K)
        n - density (cm-3)        
        z - altitude (km)        
        dz - vertical thickness (km)        
        N_int - column density, (cm-2) calculated using trapezoid rule 
        N_bar - column density, (cm-2) calculated using mean layer values 
        kmamg - column density, (km amagat) unit conversion of N_int 
    """
    
    def linear_pressure():
        levels = np.linspace(min(HASI['p']), max(HASI['p']), nlev)
        return levels

    def log_pressure():
        levels = 10**np.linspace(np.log10(min(HASI['p'])), np.log10(max(HASI['p'])), nlev)
        return levels

    def split_at_tropopause(): 
        tropopause = 300 # [mbar]
        levels = np.unique(np.append(
                    np.linspace(min(HASI['p']), tropopause, (nlev/2)+1),
                    np.linspace(tropopause, max(HASI['p']), (nlev/2)+1) ))
        return levels
                
    def equal_columns():
        dz = HASI['z'][0:-1] - HASI['z'][1:]
        N_cum = np.cumsum(dz*1e5*HASI['n'][0:-1])
        dN = N_cum[-1]/nlay
        levels = np.append(np.append(HASI['p'][0], 
                                     HASI['p'][fi(N_cum, dN*np.arange(1,nlay))]), 
                           HASI['p'][-1])
        return levels

    altitude_scheme = {'linear_pressure':linear_pressure,
                       'log_pressure':log_pressure,
                       'split_at_tropopause':split_at_tropopause,
                       'equal_columns':equal_columns,
                       }
    
    levels = altitude_scheme[method]()
        
    nlev = len(levels)
    nlay = nlev-1    
    layers = {'p':np.ndarray(nlay), 'p_min':np.ndarray(nlay), 'p_max':np.ndarray(nlay),
              'T':np.ndarray(nlay), 'T_min':np.ndarray(nlay), 'T_max':np.ndarray(nlay),
              'z':np.ndarray(nlay), 'z_min':np.ndarray(nlay), 'z_max':np.ndarray(nlay),
              'n':np.ndarray(nlay), 'n_min':np.ndarray(nlay), 'n_max':np.ndarray(nlay),
              'dz':np.ndarray(nlay),'N_bar':np.ndarray(nlay), 'N_int':np.ndarray(nlay),
              'kmamg':np.ndarray(nlay), 'method':method,
              }
    
    for i in range(nlay):
        p_min = levels[i] 
        p_max = levels[i+1] 
        
        layers['p_min'][i] = p_min
        layers['p_max'][i] = p_max
        layers['p'][i] = np.sqrt(p_min*p_max)

        layers['T_min'][i] = temperature_at_pressure(p_min)
        layers['T_max'][i] = temperature_at_pressure(p_max)
        layers['T'][i] = np.sqrt(layers['T_min'][i]*layers['T_max'][i])

        layers['n_min'][i] = density_at_pressure(p_min)
        layers['n_max'][i] = density_at_pressure(p_max)
        layers['n'][i] =  np.sqrt(layers['n_min'][i]*layers['n_max'][i])
        
        layers['z_min'][i] = altitude_at_pressure(p_max)
        layers['z_max'][i] = altitude_at_pressure(p_min)
        layers['dz'][i] = layers['z_max'][i] - layers['z_min'][i]
        layers['z'][i] = (layers['z_max'][i]+layers['z_min'][i])/2.
        
        layers['N_bar'][i] = layers['n'][i]*layers['dz'][i]*1e5

        n_sub_levels = 20
        sub_levels = np.linspace(p_min, p_max, n_sub_levels)
        z = altitude_at_pressure(sub_levels)        
        n = density_at_pressure(sub_levels)
        layers['N_int'][i] = np.trapz(n, -z*1e5)

        layers['kmamg'][i] = layers['N_int'][i]/amg/1e5
                
    model = {'nlev':nlev, 'nlay':nlay, 'levels':levels, 'layers':layers,}
    
    return model

def plot_structure(model, ax=None, logP=False):
    """Plot of HASI structure with model levels."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    axfi = lambda array, v : (np.abs(array-v)).argmin()

    levels = model['levels']
    
    if not ax: fig, ax = plt.subplots()
    ax.set_xlabel('temperature, T (K)')
    if logP:
        alt_km  = lambda p, pos : "{:4.1f}".format(HASI['z'][axfi(HASI['p'], 10**p)])
        formatter = FuncFormatter(alt_km)
        ax.plot(HASI['T'], np.log10(HASI['p']), 'k')
        ax.set_ylabel('pressure, log p (mbar)')
        ax.set_ylim(np.log10([max(HASI['p']), min(HASI['p'])]))
        for p in levels: ax.plot([min(HASI['T']),max(HASI['T'])],[np.log10(p),np.log10(p)], 'k--')
    else:
        alt_km  = lambda p, pos : "{:4.1f}".format(HASI['z'][axfi(HASI['p'], p)])
        formatter = FuncFormatter(alt_km)
        ax.plot(HASI['T'], HASI['p'], 'k')
        ax.set_ylabel('pressure, (mbar)')
        ax.set_ylim([max(HASI['p']), min(HASI['p'])])
        for p in levels: ax.plot([min(HASI['T']),max(HASI['T'])],[p,p], 'k--')
       
    axr = ax.twinx()
    axr.set_ylim(ax.get_ylim())
    axr.set_ylabel('altitude, z (km)')
    axr.yaxis.set_major_formatter(formatter)

def show_atmosphere_structures(model, filename=None):
    """Figure for comparing the layer distributions on 
    log pressure and linear pressure scales against the 
    temperature vertical profile."""

    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(19,11))
    fig.subplots_adjust(hspace=0.2, wspace=0.3)
    for i,method in enumerate(['linear_pressure', 'split_at_tropopause','log_pressure', ]):
        axs[0,i].set_title(method)
        plot_structure(set_HASI_structure(nlev=model['nlev'], method=method), 
                       logP=True, ax=axs[0,i])
        plot_structure(set_HASI_structure(nlev=model['nlev'], method=method), 
                       logP=False, ax=axs[1,i])
    if filename: savefig(filename)

def structure_table_with_boundaries(layers, latex_file=None):
    """Create latex table of layer propoerties."""    
    hdr = r"""
\begin{table}[ht]
\caption{Atmospheric Structure: Properties at Layer Boundaries and Column Comparison} 
\begin{center}
\begin{tabular}{c|ccc|ccc|ccc|ccc|cc} 
\hline\hline
"""    
    hdr += "{:^s} & {:^24s} &{:^24s} &{:^24s} &{:^24s} &{:^24s} \\\\".format(
           'Layer',
           '\multicolumn{3}{c}{Pressure (mbar)}', 
           '\multicolumn{3}{c}{Temperature (K)}', 
           '\multicolumn{3}{c}{Altitude (km)}',
           '\multicolumn{3}{c}{Density (cm$^{{-3}}$)}',
           '\multicolumn{2}{c}{Column (cm$^{{-2}}$)}',
            )
    subhdr = "{:^11s}&{:^3s}&{:^11s}&".format(r'$z_{\rm min}$', '', r'$z_{\rm max}$')
    line   = "{:3.0f} & "+r"{:^8.1f} & {{\bf {:^5.1f}}} &{:^8.1f} &"*3+ \
             "{:^9.2e}& {{\\bf {:7.2e}}} &{:^9.2e}  & {{\\bf {:^9.2e}}} &{:^9.2e}\\\\ \n"
 
    tbl_txt = hdr+'\n'+"&"+(subhdr*4)+" $N_{{\\rm int}}$ & $N_{{\\rm mean}}$ \\\\ \\hline \n"
    for i in range(len(layers['p'])): 
        tbl_txt += line.format(i+1, 
          layers['p_min'][i], layers['p'][i], layers['p_max'][i],
          layers['T_min'][i], layers['T'][i], layers['T_max'][i], 
          layers['z_min'][i], layers['z'][i], layers['z_max'][i], 
          layers['n_min'][i], layers['n'][i], layers['n_max'][i],
          layers['N_int'][i], layers['N_bar'][i],
          ) 

    tbl_txt += """\hline\n\end{tabular}\n\end{center}\end{table}"""    

    if latex_file:
        with open(latex_file, "w") as f:
            f.write(tbl_txt)
            
    return tbl_txt

def structure_table(layers, latex_file=None):
    """Create latex table of layer propoerties."""    
    hdr = r"""
\begin{table}[ht]
\caption{Atmospheric Structure}
\begin{center}
\begin{tabular}{*{7}{c}} 
\hline\hline
"""
    hdr += " {:s} & {:s} & {:s} & {:s} & {:s} & {:s} \\\\ \n".format(
           'Layer', 'Pressure', 'Altitude',
           'Temperature', 'Density', '\multicolumn{2}{c}{Column}',
            )

    subhdr = " {:s} & {:s} & {:s} & {:s} & {:s} & {:s} & {:s} \\\\ \hline \n".format(
             ' ', '(mbar)', '(km)', '(K)', '(cm$^{{-3}}$)', '(cm$^{{-2}}$)', '(km amg)',
            )

    line   = "{:3.0f}"+"&{:^8.1f}"*3+"&{:^9.2e}"*2+"&{:^8.2f}"+" \\\\ \n"
 
    tbl_txt = hdr+subhdr
    for i in range(len(layers['p'])):
        tbl_txt += line.format(i+1, 
          layers['p'][i],
          layers['z'][i], 
          layers['T'][i], 
          layers['n'][i], 
          layers['N_int'][i],
          layers['kmamg'][i],
          ) 

    tbl_txt += """\hline\n\end{tabular}\n\end{center}\end{table}"""    

    if latex_file:
        with open(latex_file, "w") as f:
            f.write(tbl_txt)
    
    return tbl_txt

def write_structure_tables(model):
    """Write latex table to output file in documentation directory."""
    structure_table(model['layers'], 
            latex_file='../doc/model_struct.tex')
    structure_table_with_boundaries(model['layers'], 
            latex_file='../doc/atm_struct_w_levels.tex')    
    return


def layer_column_density(p_min, p_max, nlev=10):
    """Sub-divide atmopheric layer into n discrete levels 
    space uniformly in log-pressure and integrate number density."""
    
    levels = 10**np.linspace(np.log10(p_max), np.log10(p_min), nlev)
    n = density_at_pressure(levels)
    z = altitude_at_pressure(levels)
    n_int = np.trapz(n,z)

    return n_int

def test_layer_column_integration():
    model = set_structure(HASI, nlev=17)
    layers = model['layers']
    km_am = layers['dz']*layers['n']*3.73e-20
    print("{:^3s}{:>10s}{:>10s}{:>8s}{:>10s}{:>10s}{:>8s}".format(
          "i", "p_min", "p_min", "N_bar", "N_int(10)", "N_int(20)","dN(%)"))
    for i in range(model['nlay']):
        km_am_int = layer_column_density(layers['p_min'][i], 
                                         layers['p_max'][i], 
                                         nlev=20)*3.73e-20
        print("{:>3d}{:10.3f}{:10.3f}{:8.3f}{:10.3f}{:10.3f}{:8.2f}".format(
              i, layers['p_min'][i], layers['p_max'][i],
              km_am[i], 
              layer_column_density(layers['p_min'][i], layers['p_max'][i])*3.73e-20,
              km_am_int,
              100*(km_am_int-km_am[i])/km_am_int, ))