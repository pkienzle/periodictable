.. _introduction:

##########################
Periodic Table of Elements
##########################

The periodictable package provides an extensible periodic table of the elements
pre-populated with data important to neutron and X-ray scattering experiments.
With the elements pacakge you can compute the scattering potential of
a compound at a given wavelength.

There is a wealth of possible information that could be stored in
such a table, and collecting it all is far beyond the scope of this project.
Instead, we provide an extensible table in which third party packages can
provide properties in addition to the base properties we define here.

.. figure:: /images/neutron_sld.png
   :width: 50%
   :alt: neutron SLD plot

   Neutron scattering potential for each element

********
Features
********

Standard properties
   Name, symbol, mass and density are built in.

Chemical Formula
   Parses chemical formulas and computes properties such as molar mass.</dd>

Isotopes
   Mass and relative abundance of isotopes are included for known isotopes.
   Formulas can include isotope composition.

Ions
   Magnetic form factors and ionic radii for various ions.

Neutron and X-ray scattering factors
   Provides neutron and wavelength dependent X-ray scattering factors for
   elements, isotopes, and formulas.

Extensible
   New properties can be added to the periodic table from outside the package.
   Specialized tables can be created with alternatives to the standard values.

Sourced
   References are available for all information in the table.
