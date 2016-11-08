atmosphere
============

A Near-IR radiative transfer model of Titan's atmosphere.

__Summary__

Jupyter notebooks describe examples of code used in the
atmosphere package. The atmospheric structure is based on observations
by Cassini Huygens. Methane opacities (k-coefficients) calculated 
for various spectral plate scales (resolutions) using
line-by-line methods and the HITRAN2012 line list. 

Description of the model published in:

	Adamkovics et al., (2016), Icarus,
	(http://dx.doi.org/10.1016/j.icarus.2015.05.023)

__Installation__

1. Check Dependencies (e.g., SWIG). Swig 3.0.7 seems to work:

    Version 3.0.7, which can be installed with brew, should work. 
    
    A previous working version is 3.0.2
   
    (caution: version 3.0.5 doesn't work)

	Example Installation:
	
		a. Download from swig.org
	
		b. ./configure --without-pcre
	
		c. sudo make

		d. sudo make install

	Version 3.0.5 creates a _cdisort.so module with "_swigconstant"
	added inconsistently. 

2. Download the source files and run:

	python setup.py install

3. Set the environment variable, RTDATAPATH

This is where the reference data files will be located.
The path environment variable be set in a Jupyter notebook, see
example notebooks. Otherwise, use something along the lines of:

      export RTDATAPATH="/path/to/your/reference/data"

The directory can be located anywhere, and will be created
if it does not exist when downloading the reference data 
in step 3.

4. Download reference data files. 

From a python interpreter or notebook:

	import atmosphere as atm
	atm.refdata.setup_directory()


__Dependencies__

Python packages used:

	numpy
	scipy
	matplotlib
	pyfits

Running in the Jupyter notebook is recommended but not necessary.
The [CDISORT](http://www.libradtran.org/bin/cdisort-2.1.3.tar.gz) 
implementation included here, provided by Tim Dowling via 
[libradtran](http://www.libradtran.org/) requires SWIG for 
compilation at install.


__Copying__

This software is licenced under the 
[GNU Generel Public License](http://www.gnu.org/licenses/gpl.txt)

__To-do__

 - parallel calculations work from notebook, but not when methods
   defined in modules. This should be cleaned up a bit. 
