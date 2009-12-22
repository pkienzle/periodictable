#!/usr/bin/env python

from setuptools import setup, find_packages
import fix_setuptools_chmod
import sys

import periodictable

if len(sys.argv) == 1:
    sys.argv.append('install')
dist = setup(
        name = 'periodictable',
        version = '0.9',
        packages = find_packages(),
        package_data = {'periodictable' : ['xsf/*.nff','read.me']},
        #data_files = periodictable.data_files(),
        install_requires = ['pyparsing'],
)

# End of file
