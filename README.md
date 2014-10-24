atmosphere
============

A Near-IR radiative transfer model of Titan's atmosphere.

__Summary__

IPython notebooks describe examples of code used in the
atmosphere package. The atmospheric structure is based on observations
by Cassini Huygens. Methane opacities (k-coefficients) calculated 
for VLT SINFONI spectral plate scale and resolution using
line-by-line methods and the HITRAN2012 line list. 

__Installation__

1. Download the source files and run:

	python setup.py install

2. Set the environment variable, RTDATAPATH

This is where the reference data files will be located.
For example, put the following line:

      export RTDATAPATH="/path/to/your/reference/data"

in the user shell configuration file, which is usually one 
of the following:

      $HOME/.bash_profile
      $HOME/.bashrc

or for cshell add the following line,

      setenv RTDATAPATH "/path/to/your/reference/data"
      
to the shell configuration file $HOME/.cshrc

The directory can be located anywhere, and will be created
if it does not exist when downloading the reference data 
in the following step.

3. Download reference data files. 

From a python interpreter or notebook:

	import atmosphere as atm
	atm.refdata.setup_directory()

__Dependencies__

1. SWIG 
2. Python packages:

	numpy
	scipy
	matplotlib
	pyfits

Running in the IPython notebook is recommended but not necessary.
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
