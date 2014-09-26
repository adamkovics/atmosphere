#!/usr/bin/env python3
"""
setup.py file for atmosphere package
"""

#from distutils.core import setup, Extension
from setuptools import setup, Extension

cdisort = Extension('_cdisort',
                     sources=['src/cdisort.i',
                              'src/cdisort.c',
                     		      'src/locate.c',
                              'src/disotest.c',
                     		  ],
                     extra_compile_args=['-fno-strict-aliasing',
                                         '-fno-common',
                                         '-dynamic',
                                         '-Qunused-arguments',
                                         ],
                    )

setup(name='atmosphere',
      version='0.0.1',
      description='Radiative transfer model for planetary atmospheres.',
      author='Mate Adamkovics',
      author_email='mate@berkeley.edu',
      url='http://astro.berkeley.edu:~/madamkov/',
      packages=['atmosphere','atmosphere/rt'],
      ext_modules=[cdisort]
      )

# INSTRUCTIONS ON PATHS TO REFERENCE DATA
#print("""
#Set the environment variable RTDATAPATH to the location where
#reference data will be stored.
#
#For example, put the following line:
#
#      export RTDATAPATH="{dir:s}"
#
#in the user bash configuration file, which is usually one 
#of the following:
#
#      $HOME/.bash_profile
#      $HOME/.bashrc
#
#or for cshell add the following line, e.g., 
#
#      setenv RTDATAPATH "{dir:s}"
#      
#to the shell configuration file $HOME/.cshrc
#
#The directory can be located anywhere, and will be created
#if it does not exist when using the following method to 
#download reference data:
#
#      atmosphere.get_reference_data()
#
#""".format(dir=datadir))
