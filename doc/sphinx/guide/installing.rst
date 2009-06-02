.. _installing:

*********
Obtaining
*********

The periodic table will be provided as an egg on PyPI, and can be obtained
simply with::

   easy_install periodictable


The source is available via svn::

   svn co svn://danse.us/common/elements/trunk periodictable
   cd periodictable
   python setup.py develop

Track updates to the package using::

   svn up

By using the 'develop' keyword on setup.py, you can modify and update
the package in place without needing to reinstall each time.

If you find you need to modify the periodic table package, please
generate a patch and send it to the danse-dev mailing list
at cacr.caltech.edu::

   svn diff > periodictable.patch

The TortoiseSVN package on Windows has similar capabilities.
