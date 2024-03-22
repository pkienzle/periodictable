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
