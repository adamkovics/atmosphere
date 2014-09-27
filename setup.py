#!/usr/bin/env python3
"""
setup.py file for atmosphere package
"""

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
                                         ],
                    )

setup(name='atmosphere',
      version='0.0.1',
      description='Radiative transfer model for planetary atmospheres.',
      install_requires=['numpy','scipy','pyfits'],
      author='Mate Adamkovics',
      author_email='mate@berkeley.edu',
      url='http://astro.berkeley.edu:~/madamkov/',
      packages=['atmosphere','atmosphere/rt'],
      ext_modules=[cdisort]
      )
