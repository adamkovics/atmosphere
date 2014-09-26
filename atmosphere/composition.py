"""
Set the gas composition of atmosphere.
"""
import numpy as np
from atmosphere.structure import pressure_at_altitude 
from atmosphere.structure import altitude_at_pressure 
from atmosphere.structure import density_at_pressure
from atmosphere.structure import amg

def set_abundances(model, trace_gas={'m_H2':0.001}):
    """Set the chemical composition in each layer for an atmosphere 
    model that is returned by the `set_HASI_structure()` method. 
    """

    model['layers'].update(
        {'m_CH4':CH4_Niemann(model['layers']['p'], pressure=True),
         'N_CH4':integrate_CH4_column(model['layers']),
         })

    set_trace_gas(model['layers'], **trace_gas)
    return model

def CH4_Niemann(z, pressure=False):
    """Return a linearly interpolated value of the CH4 mixing
    ratio from the Niemann et al., 2005 GCMS data, assuming 
    uniform mixing ratio above and below measurement limits.
    
    Input z is altitude (km) unless pressure=True, 
    in which case input z is in units of (mbar) .    

    Note: scipy version 0.14.0 (or greater?) is required,
          otherwise the bounds error returns nans.
    """

    from scipy.interpolate import interp1d

    Niemann_data = [[134,114.5,85.5,43,37.5,32.5,
                 27,23,19,16,12.5,10.5,7,4,2],             
                [.014,.014,.014,.016,.016,.015,.017,.020,
                 .022,.025,.032,.038,.046,.049,.049]]
    
    if pressure:
        f_CH4 = interp1d(pressure_at_altitude(Niemann_data[0]), 
                         Niemann_data[1], 
                         bounds_error=False)
        m_CH4 = f_CH4(z)    
    else:
        f_CH4 = interp1d(np.log10(Niemann_data[0]), 
                 Niemann_data[1], 
                 bounds_error=False)
        m_CH4 = f_CH4(np.log10(z))

    if m_CH4.size == 1 :
        if np.isnan(m_CH4) :
            if pressure:
                if (np.array(z) > 1329):
                    m_CH4 = np.array(0.049)
                    print('Input pressure exceeds surface pressure.')
                else:
                    m_CH4 = np.array(0.0014)
            else:
                if (np.array(z) < 2):
                    m_CH4 = np.array(0.049)
                else:
                    m_CH4 = np.array(0.0014)
    else :
        if pressure:
            m_CH4[np.where(np.isnan(m_CH4) & (np.array(z) > 1329) )] = 0.049
            m_CH4[np.where(np.isnan(m_CH4) & (np.array(z) < 3.74) )] = 0.014
        else:
            m_CH4[np.where(np.isnan(m_CH4) & (np.array(z) < 2) )] = 0.049
            m_CH4[np.where(np.isnan(m_CH4) & (np.array(z)>134) ) ] = 0.014

    return m_CH4

def integrate_CH4_column(layers, nlev=10, CH4_scale=None, verbose=False):
    """Interpolate the CH4 column in layers using nlev sublevels
       and return columns in units of km amagats."""
    
    nlay = len(layers['p'])
    N_CH4 = np.ndarray(nlay)
    if verbose:        
        header = '{:>2s}'+'{:^12s}'*7
        print(header.format('i','N_CH4','min(m_CH4)','max(m_CH4)'
                            ,'min(n)','max(n)','N','dz'))
    for i in range(nlay):
        p_min = layers['p_min'][i]        
        p_max = layers['p_max'][i]
        levels = np.linspace(p_min, p_max, nlev)
        z = altitude_at_pressure(levels)
        n = density_at_pressure(levels)
        m_CH4 = CH4_Niemann(z)
        if CH4_scale:
            if np.max(z) < 30:
                m_CH4 *= CH4_scale
        N_CH4[i] = np.trapz(n*m_CH4, -z*1e5)/amg/1e5
        line = "{:>2d}"+"{:^12.2e}"*7
        if verbose:
            print(line.format(i, N_CH4[i], min(m_CH4), max(m_CH4), 
                          min(n), max(n), np.trapz(n, -z*1e5), max(z)- min(z)))
    return N_CH4

def structure_table_with_composition(layers, latex_file=None):
    """Create latex table of layer propoerties."""    
    hdr = r"""
\begin{table}[ht]
\caption{Atmospheric Structure and Composition}
\begin{center}
\begin{tabular}{*{8}{c}} 
\hline\hline
"""
    hdr += " {:s} & {:s} & {:s} & {:s} & {:s} & {:s} & {:s} \\\\ \n".format(
           'Layer', 'Pressure', 'Altitude', 'Temperature', 'Density', 
           '\multicolumn{2}{c}{Total Column}', 'CH$_4$ Column',
            )

    subhdr_fmt = " {:s} &"*7+" {:s} \\\\ \hline \n"
    
    subhdr = subhdr_fmt.format(' ','(mbar)','(km)',
            '(K)','(cm$^{{-3}}$)','(cm$^{{-2}}$)','(km amg)','(km amg)',
            )

    line   = "{:3.0f}"+"&{:^8.1f}"*3+"&{:^9.2e}"*2+"&{:^8.2f}"+"&{:^8.4f}"+" \\\\ \n"
 
    tbl_txt = hdr+subhdr
    for i in range(len(layers['p'])):
        tbl_txt += line.format(i+1, 
          layers['p'][i],
          layers['z'][i], 
          layers['T'][i], 
          layers['n'][i], 
          layers['N_int'][i],
          layers['kmamg'][i],
          layers['N_CH4'][i],
          ) 

    tbl_txt += """\hline\n\end{tabular}\n\end{center}\end{table}"""    

    if latex_file:
        with open(latex_file, "w") as f:
            f.write(tbl_txt)
    
    return tbl_txt

def write_structure_tables_comp(model):
    """Write latex table to output file in documentation directory."""
    structure_table_with_composition(model['layers'], 
            latex_file='../doc/atm_struct_chem.tex')

def set_trace_gas(layers, m_H2=0.001):
    """Define the columns of H2 and N2 using an assumed 
    hydrogen mixing ratio, m_H2."""
    layers.update({'m_H2':m_H2,
                   'm_N2':1. - layers['m_CH4'] - m_H2})
    layers.update({
        'N_H2':layers['kmamg']*layers['m_H2'],
        'N_N2':layers['kmamg']-layers['N_CH4']-layers['m_H2']*layers['kmamg'],
        })


def show_methane_profile(model, ax=None):
    """Plot methane mixing ratio with atmosphere levels
    and interpolation."""

    from atmosphere.structure import HASI

    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    axfi = lambda array, v : (np.abs(array-v)).argmin()
    levels = model['levels']

    Niemann_data = [[134,114.5,85.5,43,37.5,32.5,
                     27,23,19,16,12.5,10.5,7,4,2],             
                    [.014,.014,.014,.016,.016,.015,.017,.020,
                     .022,.025,.032,.038,.046,.049,.049]]
    if not ax: 
        fig, axs = plt.subplots(nrows=2, figsize=(5,7),)
        fig.subplots_adjust(hspace=0)

    z = np.logspace(0,np.log10(1470),200)
    ax = axs[0]
    ax.plot(100*np.array(Niemann_data[1]), 
            pressure_at_altitude(Niemann_data[0]), 'ko')
    for p in model['levels']: ax.plot( ax.get_xlim(), [p,p], 'k--')
    ax.plot(100*CH4_Niemann(z, pressure=True), z, 'k-x' )
    ax.set_ylim(300,2.75)
    ax.get_xaxis().tick_top()  
    ax.get_xaxis().set_ticks([])
    ax.set_ylabel('stratosphere, pressure (mbar)')
    axr = ax.twinx()
    axr.set_ylim(ax.get_ylim())
    axr.set_ylabel('altitude, z (km)')
    
    alt_km  = lambda p, pos : "{:4.1f}".format(HASI['z'][axfi(HASI['p'], p)])
    formatter = FuncFormatter(alt_km)
    axr.yaxis.set_major_formatter(formatter)
        
    ax = axs[1]
    ax.plot(100*np.array(Niemann_data[1]), 
            pressure_at_altitude(Niemann_data[0]), 'ko')
    for p in model['levels']: ax.plot( ax.get_xlim(), [p,p], 'k--')
    ax.plot(100*CH4_Niemann(z, pressure=True), z, 'k-x' )

    axr = ax.twinx()
    ax.set_ylim(1466,300) ; axr.set_ylim(1466,300)
    axr.set_ylabel('altitude, z (km)')
    formatter = FuncFormatter(alt_km)
    axr.yaxis.set_major_formatter(formatter)
    
    ax.get_xaxis().tick_bottom()  
    ax.set_ylabel('troposphere, pressure (mbar)')
    ax.set_xlabel('methane mixing ratio (%)')
    
def show_temperature_profile(model, ax=None):
    """Plot methane mixing ratio with atmosphere levels
    and interpolation."""

    from matplotlib.ticker import FuncFormatter
    axfi = lambda array, v : (np.abs(array-v)).argmin()
    levels = model['levels']
    
    if not ax: 
        fig, axs = subplots(nrows=2, figsize=(5,7),)
        fig.subplots_adjust(hspace=0)

    z = logspace(0,np.log10(1470),200)
    ax = axs[0]
    ax.plot(model['layers']['T'], model['layers']['p'], 'ko')
    xlim = (60,140)
    for p in model['levels']: ax.plot( xlim, [p,p], 'k--')
    ax.plot(HASI['T'], HASI['p'],'k-' )
    ax.set_ylim(300,2.75)
    ax.get_xaxis().tick_top()  
    ax.get_xaxis().set_ticks([])
    ax.set_ylabel('stratosphere, pressure (mbar)')
    ax.set_xlim(xlim)
    axr = ax.twinx()
    axr.set_ylim(ax.get_ylim())
    axr.set_ylabel('altitude, z (km)')
    
    alt_km  = lambda p, pos : "{:4.1f}".format(HASI['z'][axfi(HASI['p'], p)])
    formatter = FuncFormatter(alt_km)
    axr.yaxis.set_major_formatter(formatter)

    ax = axs[1]
    ax.plot(model['layers']['T'], model['layers']['p'], 'ko')
    for p in model['levels']: ax.plot( xlim, [p,p], 'k--')
    ax.plot(HASI['T'], HASI['p'],'k-' )
    ax.set_ylim(1466,300)
    ax.set_xlim(xlim)
    ax.get_xaxis().tick_bottom()  
    ax.set_ylabel('troposphere, pressure (mbar)')
    ax.set_xlabel('temperature (K)')

    axr = ax.twinx()
    ax.set_ylim(1466,300) ; axr.set_ylim(1466,300)
    axr.set_ylabel('altitude, z (km)')
    formatter = FuncFormatter(alt_km)
    axr.yaxis.set_major_formatter(formatter)