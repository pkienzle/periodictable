#!/usr/bin/env python

"""
Usage:

./test.py
    - run all tests

./test.py --with-coverage
    - run all tests with coverage report
"""

import os, sys
import nose

# Check that we are running from the root.
path = os.getcwd()
assert os.path.exists(os.path.join(path,'periodictable','nsf.py'))

# Make sure that we have a private version of mplconfig
mplconfig = os.path.join(os.getcwd(),'.mplconfig')
os.environ['MPLCONFIGDIR'] = mplconfig
os.putenv('MPLCONFIGDIR',mplconfig)
if not os.path.exists(mplconfig): os.mkdir(mplconfig)

# Run the source tests with periodictable on the path.  By manipulating
# the path in this way, we can test without having to build and install.
# We are adding doc/sphinx to the path because the periodic table extension
# doctests need to be able to find the example extensions.
sys.path.insert(0,path)
sys.path.insert(1,os.path.join(path,'doc','sphinx'))
nose_args = [__file__,'-v','--with-doctest','--doctest-extension=.rst',
             '--cover-package=periodictable']
nose_args += sys.argv[1:]  # allow coverage arguments
nose_args += ['test','periodictable','doc/sphinx/guide']
if not nose.run(argv=nose_args): sys.exit(1)

# Run isolated tests in their own environment.  In this case we will have
# to set the PYTHONPATH environment variable before running since it is
# happening in a separate process.
if 'PYTHONPATH' in os.environ:
    PYTHONPATH = path + ":" + os.environ['PYTHONPATH']
else:
    PYTHONPATH = path
os.putenv('PYTHONPATH', PYTHONPATH)
for p in ("nsfd2o", ):
    ret = os.system(" ".join( (sys.executable, "test/test_%s.py"%p) ))
    if ret != 0: sys.exit()
