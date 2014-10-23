"""
Set the gas composition of atmosphere.
"""
import numpy as np
from atmosphere.structure import pressure_at_altitude 
from atmosphere.structure import altitude_at_pressure 
from atmosphere.structure import density_at_pressure
from atmosphere.structure import amg

import io

data = """
# time   CH4    sigma  alt     P        T 
#
# Methane Mole Fraction from Huygens GCMS
#   Table 1 of Niemann et al., 2010, JGR
#   http://dx.doi.org/10.1029/2010JE003659
#
#  Units for columns are:
#    Time from entry, time, (seconds)
#    CH4 mole fraction, CH4, (%)
#    Standard Deviation in CH4, sigma, (%)
#    Altitude from surface, alt (km)
#    Pressure, P (hPa)
#    Temperature, T (K)
#
 193   1.49   1.32   139.8      3.3   162.8
 285   1.53   0.79   135.5      3.7   161.1
 344   1.49   1.09   132.8      4.0   161.0
 443   1.51   1.00   128.6      4.4   157.2
 565   1.51   1.32   123.7      5.0   156.1
 693   1.49   0.94   118.9      5.7   155.0
 813   1.46   1.01   114.7      6.4   152.6
 905   1.48   0.57   111.2      6.9   148.7
 965   1.46   0.89   106.5      7.6   148.4
1061   1.45   0.80    99.5      9.6   144.8
1183   1.46   0.74    91.8     12.2   139.6
1310   1.46   0.72    84.8     15.3   132.7
1425   1.47   0.60    79.2     18.5   125.5
1510   1.49   0.79    75.5     21.2   118.9
2846   1.50   1.83    44.6    112.8    70.5
3226   1.53   2.64    39.6    155.2    70.7
3683   1.64   2.49    34.4    214.4    71.4
4259   1.80   2.10    28.9    302.4    73.0
4766   2.02   2.32    24.6    390.7    74.9
5210   2.33   1.69    21.2    477.0    76.8
5617   2.60   1.68    18.3    563.2    78.5
5874   2.83   1.62    16.6    620.9    79.6
6065   3.02   1.82    15.3    665.8    80.4
6310   3.36   1.96    13.8    724.9    81.4
6715   3.97   2.82    11.3    827.4    83.2
7242   4.99   2.56     8.3    969.9    85.4
7527   5.36   1.48     6.7   1051.4    86.8
7756   5.55   2.01     5.5   1118.4    87.8
7956   5.75   2.20     4.5   1178.3    88.8
8214   5.76   1.89     3.2   1257.2    90.0
8457   5.72   2.16     2.0   1333.7    91.2
8644   5.73   1.48     1.1   1393.8    92.3
8807   5.67   1.50     0.3   1446.7    93.1
"""

Niemann = np.genfromtxt(io.BytesIO(data.encode()), names=True)

def set_abundances(model, trace_gas={'m_H2':0.001}):
    """Set the chemical composition in each layer for an atmosphere 
    model that is returned by the `set_HASI_structure()` method. 
    """

    model['layers'].update(
        {'m_CH4':CH4_Niemann(model['layers']['z']),
         'N_CH4':integrate_CH4_column(model['layers']),
         })

    set_trace_gas(model['layers'], **trace_gas)
    return model

def CH4_Niemann(z, pressure=False):
    """Return a linearly interpolated value of the CH4 mole fraction
    from Huygens GCMS measurement. Niemann et al., 2010.
    
    Input z is altitude (km) unless pressure=True, 
    in which case input z is in units of (mbar) .    
    """

    if pressure:
        m_CH4 = np.interp(z, Niemann['P'], Niemann['CH4']/100.)
    else:
        m_CH4 = np.interp(np.log10(z), 
                          np.log10(Niemann['alt']), 
                          Niemann['CH4']/100.,)

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


def show_methane_profile(model, ax=None, old_data=False):
    """Plot methane mixing ratio with atmosphere levels
    and interpolation."""

    from atmosphere.structure import HASI

    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    axfi = lambda array, v : (np.abs(array-v)).argmin()
    levels = model['levels']

    Niemann_2005 = [[134,114.5,85.5,43,37.5,32.5,
                     27,23,19,16,12.5,10.5,7,4,2],             
                    [.014,.014,.014,.016,.016,.015,.017,.020,
                     .022,.025,.032,.038,.046,.049,.049]]
    if not ax: 
        fig, axs = plt.subplots(nrows=2, figsize=(5,7),)
        fig.subplots_adjust(hspace=0)

    z = np.logspace(0,np.log10(1470),200)
    ax = axs[0]
    ax.plot(100*CH4_Niemann(z, pressure=True), z, 'k-x', label='interp.')
    ax.plot(Niemann['CH4'], Niemann['P'], 'ko', label='Niemnan 2010')
    if old_data:
        ax.plot(100*np.array(Niemann_2005[1]), 
                pressure_at_altitude(Niemann_2005[0]), 
                'ko', alpha=0.3, label='Nieman 2005')
    for p in model['levels']: ax.plot( ax.get_xlim(), [p,p], 'k--')
    ax.set_ylim(300,2.75)
    ax.get_xaxis().tick_top()  
    ax.get_xaxis().set_ticks([])
    ax.set_ylabel('stratosphere, pressure (mbar)')
    ax.legend(loc=1)
    axr = ax.twinx()
    axr.set_ylim(ax.get_ylim())
    axr.set_ylabel('altitude, z (km)')
    
    alt_km  = lambda p, pos : "{:4.1f}".format(HASI['z'][axfi(HASI['p'], p)])
    formatter = FuncFormatter(alt_km)
    axr.yaxis.set_major_formatter(formatter)
        
    ax = axs[1]
    ax.plot(100*CH4_Niemann(z, pressure=True), z, 'k-x' )
    ax.plot(Niemann['CH4'], Niemann['P'], 'ko', label='Niemnan 2010')
    if old_data:
        ax.plot(100*np.array(Niemann_2005[1]), 
                pressure_at_altitude(Niemann_2005[0]), 'ko', alpha=0.3)
    for p in model['levels']: ax.plot( ax.get_xlim(), [p,p], 'k--')

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