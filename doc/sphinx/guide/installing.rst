.. _installing:

**********
Installing
**********

This tutorial will walk you through the process of installing Periodic Table. 
To follow, you only need two basic things:

* A working Python installation, version 2.4 or higher.
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

If you find you need to modify the periodic table package, please
generate a patch and send it to the `DANSE <http://danse.us>`_ 
Project mailing list at danse-dev@cacr.caltech.edu::

   svn diff > periodictable.patch

Windows user can use `TortoiseSVN <http://tortoisesvn.tigris.org/>`_ 
package which provides similar operations.
