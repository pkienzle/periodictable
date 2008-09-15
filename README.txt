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

