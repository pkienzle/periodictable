#!/usr/bin/env python

from setuptools import setup, find_packages
import fix_setuptools_chmod
import sys

import periodictable

if len(sys.argv) == 1:
    sys.argv.append('install')
dist = setup(
        name = 'periodictable',
        version = periodictable.__version__,
        author='Paul Kienzle',
        author_email='pkienzle@gmail.com',
        url='http://www.reflectometry.org/danse/elements.html',
        description='Extensible periodic table of the elements',
        long_description=open('README.txt').read(),
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: Public Domain',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Physics',
            ],
        packages = find_packages(),
        package_data = {
            'periodictable' :
                ['xsf/*.nff', 'xsf/f0_WaasKirf.dat', 'xsf/read.me'],
        },
        #data_files = periodictable.data_files(),
        install_requires = ['pyparsing'],
)

# End of file
