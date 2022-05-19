=========================================
Extensible periodic table of the elements
=========================================

This package provides a periodic table of the elements with
support for mass, density and xray/neutron scattering information.

Masses, densities and natural abundances come from the
NIST Physics Laboratory, but do not represent a critical
evaluation by NIST scientists.

Neutron scattering calculations use values collected by the
Atomic Institute of the Austrian Universities.  These values
do corresponding to those from other packages, though there
are some differences depending to the tables used.  Bound
coherent neutron scattering for gold in particular is significantly
different from older value: 7.63(6) as measured in 1974
compared to 7.90(7) as measured in 1990.

X-ray scattering calculations use a combination of empirical and
theoretical values from the LBL Center for X-ray Optics.  These
values differ from those given in other sources such as the
International Tables for Crystallography, Volume C, and so may
give different results from other packages.

Install using::

    pip install periodictable

There are several web interfaces to this package:

* `NIST Center for Neutron Research (NCNR) <https://www.ncnr.nist.gov/resources/activation>`_
* `Forschungs-Neutronenquelle Heinz Maier-Leibnitz (FRM II) <https://webapps.frm2.tum.de/intranet/activation/>`_
* `SLD Calulator <https://sld-calculator.appspot.com/>`_

Documentation is `available online <https://periodictable.readthedocs.io>`_.

Source links:

* https://github.com/pkienzle/periodictable
* https://github.com/scattering/activation

|CI| |RTD| |DOI|

.. |CI| image:: https://github.com/pkienzle/periodictable/workflows/Test/badge.svg
   :alt: Build status
   :target: https://github.com/pkienzle/periodictable/actions

.. |DOI| image:: https://zenodo.org/badge/1146700.svg
   :alt: DOI tag
   :target: https://zenodo.org/badge/latestdoi/1146700

.. |RTD| image:: https://readthedocs.org/projects/periodictable/badge/?version=latest
   :alt: Documentation status
   :target: https://periodictable.readthedocs.io/en/latest/?badge=latest

Known issues
============

* Incoherent scattering computed for contrast matched mixture in D2O_sld(),
  differs from the value that would be computed for a compound with the same
  isotope proportions and density computed in neutron_sld(). This may change
  in a future release.

* The mass and composition tables are out of date. This package uses tables
  from 1997 but IUPAC produced new tables in 2009.

* Incoherent scattering calculations for energy-dependent rare earth elements
  is underestimated. The calculation requires bound incoherent scattering
  lengths (b_i) and the bound coherent lengths (b_c), but only b_c is
  included.

Change history
==============

1.6.1 2022-05-18
----------------

Modified:

* Calculate decay time correctly in the presence of significant long-lived
  activation.
* Calculate b mode activation correctly for Eu[151] => Gd[152]

1.6.0 2021-04-21
----------------

New:

* Add energy dependence for rare earths (Lynn and Seeger, 1990).

Modified:

* Use complex b_c when computing the coherent cross section, leading to
  correct values of sigma_c and sigma_i for materials with large absorption.
  With this change the tabulated values for B[10] are now shown to be
  self-consistent within a few percent.

Breaking changes:

* Neutron scattering factors are returned with one value for each wavelength
  even for energy independent elements. Previous versions returned a scalar
  if the returned value was identical for each wavelength.

1.5.3 2020-11-04
----------------

Breaking changes:

* Fix calculation of contrast match points for biomolecules. The old
  formula used the density of H2O for the D2O sld calculation.
* Modify biomolecule support to use H[1] rather than T for labile hydrogen.
  This will result in less error when the labile formula is used in lieu
  of the natural formula or the contrast-matched formula, and make it more
  obvious from glancing at the formula that labile hydrogen is present.
* Modify *fasta.Molecule* attributes, dropping *Hmass* and *Hsld*. *Hnatural*
  has been moved to *natural_formula*. The formula with labile hydrogen is
  stored in *labile_formula*, as well as *formula* as before.

New:

* Add *replace()* method to formula to allow isotope substitution.
* Add *nsf.D2O_match()* and *nsf.D2O_sld()* functions.

Modified:

* Neutron wavelength now defaults to 1.798 A when wavelength and energy are
  both None in *neutron_sld()* and *neutron_scattering()* rather than
  throwing an assertion error.
* *table* can be passed to neutron sld calculators as the source of isotope
  information when parsing the chemical formula.
* Switch unit test framework from nose to pytest.
* Update docs.

1.5.2 2019-11-19
----------------

Modified:

* Carbon density changed from 2.1 to 2.2 to match CXRO, CRC and RSC. The NIST
  X-ray attenuation tables use 2.26; the Handbook of Mineralogy has 2.09-2.23.
  The Neutron Data Booklet gave the value as 1.9-2.3, and 2.1 was chosen
  from this range.  The remaining density will continue to use values from the
  Neutron Data Booklet, which cites CRC as the primary source.
* Updated references.

1.5.1 2019-09-09
----------------

Modified:

* fasta uses natural abundance of H for biomolecule when computing the
  D2O contrast match rather than the biomolecule with pure H[1].
* remove half-life units from column header in activation table since
  each row gives its own units.

1.5.0 2017-05-11
----------------

New:

* mixture by mass and volume, e.g., 5 g NaCl // 50 mL H2O@1
* multilayer materials, e.g., 5 um Si // 3 nm Cr // 8 nm Au
* add support for bio molecules with labile hydrogens
* update list of possible oxidation states to include rare states

Modified:

* fixed computation of incoherent cross section so it is consistent with
  coherent cross section and total cross section


1.4.1 2014-02-04
----------------

Modified:

* default density is now the isotopic density rather than the natural density

1.4.0 2013-12-20
----------------

* support python 3.3

1.3.10 2013-10-25
-----------------

Modified:

* fix activation calculation to ignore fast neutrons in thermal environment
* add emission spectra for remaining elements above neon

1.3.9 2013-04-23
----------------

Modified:

* Update requirements to pyparsing<2.0.0 (we don't support python 3 yet)

1.3.8 2013-04-08
----------------

New:

* formula parser supports density spec and mix by weight/mix by volume

Modified:

* py2exe/py2app wrapping now includes missing activation.dat
* skipping bad 1.3.7 build which didn't include all changes

1.3.6 2013-03-05
----------------

New:

* add activation decay time to neutron activation calculator

Modified:

* Change neutron scattering calculations for incoherent cross section
  to be the linear combination of the incoherent cross sections of the
  individual atoms rather than total cross section minus the coherent
  cross section.  Penetration depth of the unscattered beam still uses
  the total cross section plus the absorption cross section.

1.3.5 2013-02-26
----------------

New:

* formulas now report charge and mass_fraction
* formula parser accepts ions as Yy{#+} or Yy[#]{#+} for isotopes
* support neutron activation calculations
* support xray refraction index and mirror reflectivity

Modified:

* update X-ray scattering tables for Zr
* adjust ion mass for number of electrons
* ions now display as Yy{#+} rather than Yy^{#+}
* fix formula.natural_density
* fix formula.hill so C,H come first
* fix element.interatomic_distance
* formula(value=...) -> formula(compound=...)

1.3  2010-12-05
---------------

New:

* mix_by_weight and mix_by_volume formula constructors
* use natural density to set density for isotope specific formulas
* add neutron_scattering function which returns xs, sld and penetration depth

Modified:

* need wavelength= or energy= for xray/neutron sld
* improved docs and testing

1.2  2010-04-28
---------------

New:

* support pickle: id(H) == id(loads(dumps(H)))
* support ions, with magnetic form factors and x-ray f0 scattering factor
* support py2exe wrappers
* allow density to be calculated from structure (bcc, fcc, hcp, cubic, diamond)
* estimate molecular volume
* support private tables with some values replaced by application

Modified:

* rename package periodictable
* rename table to periodictable.elements
* neutron sld returns real and imaginary coherent and incoherent
  instead of coherent, absorption and incoherent
* bug fix: sld for H[2] was wrong when queried before sld for H.
* remove CrysFML ionic radius definitions

1.1  2009-01-20
---------------

Modified:

* Restructure package, separating tests into different directory
* When defining table extensions, you should now do::

      from elements.core import periodic_table, Element, Isotope

  rather than::

      from elements import periodic_table
      from elements.elements import Element, Isotope

