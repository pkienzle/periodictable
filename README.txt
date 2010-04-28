Extensible periodic table of the elements

Release 1.2, April 28, 2010

This package provides a periodic table of the elements with
support for mass, density and xray/neutron scattering information.

Masses, densities and natural abundances come from the
NIST Physics Laboratory, but do not represent a critical
evaluation by NIST scientists.

Neutron scattering calculations use values collected by the 
Atomic Institute of the Austrian Universities.  These values
do corresponding to those from other packages, though there 
are some differences according to the tables used.  Bound 
coherent neutron scattering for gold in particular is significantly 
different from the old value: 7.63(6) as measured in 1974
compared to 7.90(7) as measured in 1990.

X-ray scattering calculations use a combination of empirical and
theoretical values from the LBL Center for X-ray Optics.  These
values differ from those given in other sources such as the
International Tables for crystallography, Volume C, and so may
give different results from other packages.


Features under consideration:

   composition and density table for common materials
   alloy density calculator
   volume fraction calculator
   density from lattice parameters
   uncertainty in sld as a function of incident energy (for X-rays),
      composition (for alloys), volume fraction (for solutions) or
      density (for pure materials).  I need this to set ranges on 
      fitting parameters.
   magnetic scattering parameters

Change history:

== 1.2  2010-04-28 ==

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

== 1.1  2009-01-20 ==

Restructure package, separating tests into different directory.

When defining table extensions, you should now do:

  from elements.core import periodic_table, Element, Isotope

rather than:

  from elements import periodic_table
  from elements.elements import Element, Isotope

