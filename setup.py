#!/usr/bin/env python

from setuptools import setup, find_packages
import fix_setuptools_chmod

dist = setup(
        name = 'periodictable',
        version = '0.9',
        packages = find_packages(),
        package_data = {'periodictable' : ['xsf/*']},
        install_requires = ['pyparsing'],
)

# End of file
