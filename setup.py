#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup, find_packages

setup(name='slitless',
      version='1.0.0',
      description='Tools for slitless spectroscopy',
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/dnidever/slitless',
      packages=['slitless'],
      package_dir={'':'python'},
      #scripts=['bin/slitless'],
      install_requires=['numpy','astropy(>=4.0)','scipy','dlnpyutils(>=1.0.3)']
      #include_package_data=True,
)
