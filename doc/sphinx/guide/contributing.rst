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
    python2.5 tests.py
    python2.6 tests.py

When all the tests run, generate a patch and send it to the 
`DANSE <http://danse.us>`_ Project mailing list at danse-dev@cacr.caltech.edu::

    git diff > periodictable.patch

Alternatively, create a project fork at github and we can pull your
changes directly from your repository.

Windows user can use `TortoiseGit <http://code.google.com/p/tortoisegit/>`_ 
package which provides similar operations.

Building the package documentation requires a working sphinx installation 
and in addition, a copy of `MathJax <http://www.mathjax.org/>`_ to view 
the equations.  Download and unzip the MathJax package into the doc/sphinx
directory to install MathJax.  You can then build the documentation as follows::

    (cd doc/sphinx && make clean html pdf)

You can see the result by pointing your browser to::

    periodictable/doc/sphinx/_build/html/index.html
    periodictable/doc/sphinx/_build/latex/PeriodicTable.pdf

As of this writing, the \\AA LaTeX command for the Angstrom symbol is not
available in the MathJax distribution. We patched jax/input/TeX/jax.js
with the additional symbol AA using::

    // Ord symbols
    S:            '00A7',
  + AA:           '212B',
    aleph:        ['2135',{mathvariant: MML.VARIANT.NORMAL}],

If you are using unusual math characters, you may need similar patches 
for your own documentation.

ReStructured text format does not have a nice syntax for superscripts and
subscripts.  Units such as |g/cm^3| are entered using macros such as
\|g/cm^3| to hide the details.  The complete list of macros is available in

        doc/sphinx/rst_prolog

In addition to macros for units, we also define cdot, angstrom and degrees 
unicode characters here.  The corresponding latex symbols are defined in 
doc/sphinx/conf.py.
