.. _introduction:

##########################
Periodic Table of Elements
##########################

The periodictable package provides an extensible periodic table with various 
properties of the elements like name, symbol, mass, desnity etc and also 
provides data important to neutron and X-ray scattering experiments. With the 
elements package you can compute the scattering potential of a compound at a
given wavelength.

There is a wealth of possible information that could be stored in
such a table, and collecting it all is far beyond the scope of this project.
Instead, we provide an extensible table in which third party packages can
provide properties in addition to the base properties we define here.

.. plot:: plots/sld_plot.py


Neutron SLD as a function of element.

********
Features
********

Standard properties
   Name, symbol, :mod:`mass <periodictable.mass>` and 
   :mod:`density <periodictable.density>` of elements are built in.

:class:`Chemical Formula <periodictable.formulas.Formula>`
   Parses chemical formula and computes properties such as molar mass.

:class:`Isotopes <periodictable.core.Isotope>`
   Mass and relative abundance of isotopes are included for known isotopes.
   Formulas can include isotope composition.

:class:`Ions <periodictable.core.Ion>`
   Magnetic form factors and ionic radii for various ions.

:mod:`Neutron <periodictable.nsf>` and :mod:`X-ray Scattering Factors <periodictable.xsf>`
   Provides neutron and wavelength dependent X-ray scattering factors for
   elements, isotopes, and formulas.

:ref:`Extensible <extending>`
   New properties can be added to the :class:`Periodic Table <periodictable.core.PeriodicTable>` 
   from outside the package.  Specialized tables can be created with 
   alternatives to the standard values.

:ref:`Data Sources <data-sources>`
   References are available for all information in the table.


