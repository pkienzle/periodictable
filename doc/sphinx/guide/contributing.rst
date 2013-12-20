.. _contributing:

********************
Contributing Changes
********************

The best way to contribute to the periodic table package is to work
from a copy of the source tree in the revision control system.

The source is available via git::

    git clone https://github.com/pkienzle/periodictable.git
    cd periodictable
    python setup.py develop

By using the *develop* keyword on setup.py, changes to the files in the
package are immediately available without the need to run setup.py
install each time.

Track updates to the original package using::

    git pull

If you find you need to modify the periodic table package, please update
the documentation and add tests for your changes.  We use doctests on all
of our examples that we know our documentation is correct.  More thorough
tests are found in test directory.  Using the the nose test package, you 
can run both sets of tests::

    easy_install nose
    python2.6 tests.py
    python2.7 tests.py
    python3.3 tests.py

When all the tests run, generate a patch and send it to the 
`DANSE <http://danse.us>`_ Project mailing list at danse-dev@cacr.caltech.edu::

    git diff > patch

Alternatively, create a project fork at github and we can pull your
changes directly from your repository.

Windows user can use `TortoiseGit <http://code.google.com/p/tortoisegit/>`_ 
package which provides similar operations.

You can then build the documentation as follows::

    easy_install sphinx
    (cd doc/sphinx && make clean html pdf)

You can see the result by pointing your browser to::

    periodictable/doc/sphinx/_build/html/index.html
    periodictable/doc/sphinx/_build/latex/PeriodicTable.pdf

ReStructured text format does not have a nice syntax for superscripts and
subscripts.  Units such as |g/cm^3| are entered using macros such as
\|g/cm^3| to hide the details.  The complete list of macros is available in

        doc/sphinx/rst_prolog

In addition to macros for units, we also define cdot, angstrom and degrees 
unicode characters here.  The corresponding latex symbols are defined in 
doc/sphinx/conf.py.
