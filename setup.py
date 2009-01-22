#!/usr/bin/env python

from setuptools import setup, find_packages
import fix_setuptools_chmod

dist = setup(
        name = 'elements',
        version = '1.1',
        packages = find_packages(),
        package_data = {'elements' : ['xsf/*']},
        install_requires = ['pyparsing'],
)

# End of file
