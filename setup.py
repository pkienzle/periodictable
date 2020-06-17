#!/usr/bin/env python
import sys
import os

from setuptools import setup, find_packages

ROOT = os.path.join(os.path.dirname(__file__))

VERSION = None
for line in open(os.path.join(ROOT, "periodictable", "__init__.py")):
    if "__version__" in line:
        VERSION = line.split('"')[1]

if len(sys.argv) == 1:
    sys.argv.append('install')

def readtext(path):
    with open(os.path.join(ROOT, path)) as fid:
        return fid.read()

setup(
    name='periodictable',
    version=VERSION,
    author='Paul Kienzle',
    author_email='pkienzle@gmail.com',
    license='public domain',
    url='https://github.com/pkienzle/periodictable',
    description='Extensible periodic table of the elements',
    long_description=readtext('README.rst'),
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
    tests_require=['pytest', 'pytest-cov'],
)

# End of file
