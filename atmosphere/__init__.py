"""

A model of Titan's atmosphere with radiative transfer solvers.

Modules:

	refdata - obtain reference data for gas and aerosols 
	structure - setup vertical structure and layers
	composition - set chemical abundances in each layer
	gas_opacity - calculate gas opacity using k-coefficients
	aerosol - determine aerosol opacity and phase functions

Sub-packages:

	rt - radiative transfer solvers including DISORT and twostream

Example:

    model = atm.structure.set_HASI_structure(
    		   	nlev=21, method='split_at_tropopause'
    		   	)
    atm.composition.set_abundances(
    	model, trace_gas={'m_H2':0.001}
    	)
    atm.gas_opacity.set_methane(
        model, kc_file=os.path.join(
        				   os.getenv('RTDATAPATH'),
                           'gas_opacity/kc_CH4.SINFONI.v08.dnu_3.0.fits')
						   )
    atm.gas_opacity.set_cia(model)
    atm.aerosol.set_opacity(model)
    DISR = atm.aerosol.fit_DISR_phases()
    atm.aerosol.set_aerosol_phase_moments(model, DISR, nmom=nmom)

Source availabe from:
	
	https://github.com/adamkovics/atmosphere

"""

try:
  from atmosphere import composition
  from atmosphere import structure
  from atmosphere import aerosol
  from atmosphere import gas_opacity
  from atmosphere import refdata
except:
  None


__all__ = ['composition',
           'structure',
           'aerosol',
           'gas_opacity',
           'refdata'
           ]

__version__ = '0.0.4'