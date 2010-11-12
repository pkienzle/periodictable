.. _installing:

**********
Installing
**********

This tutorial will walk you through the process of installing Periodic Table. 
To follow, you need some basic things:

* A working Python installation (tested on 2.5 and 2.6)
* A working numpy package
* Easy_install module, if you don't have easy_install installed on your 
  system, download `here <http://pypi.python.org/pypi/setuptools#files>`_.

The periodic table will be provided as an egg on `PyPI <http://pypi.python.org/pypi>`_, 
and can be obtained simply with::

    [coming soon]  easy_install periodictable

The source is available via svn::

    svn co svn://danse.us/common/elements/trunk periodictable
    cd periodictable
    python setup.py develop

Track updates to the package using::

    svn update

By using the *develop* keyword on setup.py, you can modify and update
the package in place without needing to reinstall each time.

If you find you need to modify the periodic table package, please update
the documenation and add tests for your changes.  We use a combination
of doctests so that we know our documented examples work as expected and
more thorough tests in tests directory.  Using the nosetests package you 
can run the test suite::

    python2.5 tests.py   
    python2.6 tests.py

When all the tests run, generate a patch and send it to the 
`DANSE <http://danse.us>`_ Project mailing list at danse-dev@cacr.caltech.edu::

    svn diff > periodictable.patch

Windows user can use `TortoiseSVN <http://tortoisesvn.tigris.org/>`_ 
package which provides similar operations.

Building the package documentation requires a working sphinx installation
and a copy of MathJax.  See the README file in the doc/sphinx directory
for details.  Once you have sphinx and MathJax set up you can build the
documentation as follows::

    (cd doc/sphinx && make clean html pdf)


