#!/usr/bin/env python
import sys
import os

from setuptools import setup, find_packages
#import fix_setuptools_chmod

version = None
for line in open(os.path.join("periodictable", "__init__.py")):
    if "__version__" in line:
        version = line.split('"')[1]

if len(sys.argv) == 1:
    sys.argv.append('install')
setup(
    name='periodictable',
    version=version,
    author='Paul Kienzle',
    author_email='pkienzle@gmail.com',
    license='public domain',
    url='http://www.reflectometry.org/danse/elements.html',
    description='Extensible periodic table of the elements',
    long_description=open('README.rst').read() if os.path.exists('README.rst') else None,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: Public Domain',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        ],
    packages=find_packages(),
    include_package_data=True,
    package_data={
        # NOTE: be sure to include files in MANIFEST.in as well
        'periodictable' :
            ['activate.dat', 'xsf/*.nff', 'xsf/f0_WaasKirf.dat', 'xsf/read.me'],
    },
    #data_files = periodictable.data_files(),
    install_requires=['pyparsing', 'numpy'],
)

# End of file
