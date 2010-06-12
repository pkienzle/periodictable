.. _formula:


********************
Chemical Composition
********************

Some properties are available for groups of elements.  Groups are specified
as a chemical formula string and either density or cell volume for the crystal
structure.  While it does not provide any information about molecular 
structure, a formula does all complete control over chemical composition. 

* Individual atoms are represented by periodic table symbol.  These are
  case sensitive, so "CO" is different from "Co".
  
* Formula strings consist of counts and atoms such as "CaCO3+6H2O".

* Groups can be separated by '+' or space, so "CaCO3 + 6H2O" works as well. 

* Groups can be defined using parentheses, such as "CaCO3 (H2O)6". 

* Parentheses can nest: "(CaCO3(H2O)6)1". 

* Isotopes are represented by index, e.g., "CaCO[18]3+6H2O". 

* Counts can be integer or decimal, e.g. "CaCO3+(3HO1.5)2".

A formula string is translated into a formula using 
:func:`periodictable.formula`.  Once the formula has been formed,
you can perform algebra on the entire formula, such as adding
together two formulas to make a more complex compound.

The following is an example of hydrated quartz:

.. doctest::

    >>> SiO2 = periodictable.formula('SiO2')
    >>> hydrated = SiO2 + periodictable.formula('3H2O')
    >>> print hydrated,'mass',hydrated.mass
    SiO2(H2O)3 mass 114.13014
    >>> rho,mu,inc = periodictable.neutron_sld('SiO2+3H2O',density=1.5,wavelength=4.75)
    >>> print hydrated,'neutron sld','%.3g'%rho
    SiO2(H2O)3 neutron sld 0.849
    >>> rho,mu = periodictable.xray_sld(hydrated,density=1.5,wavelength=Cu.K_alpha)
    >>> print hydrated,'X-ray sld','%.3g'%rho
    SiO2(H2O)3 X-ray sld 13.5

Formulas are parsed from strings using the following grammar::

        number    :: [1-9][0-9]*
        fraction  :: ( | '0' | number) '.' [0-9]*
        count     :: number | fraction | ''
        symbol    :: [A-Z][a-z]*
        isotope   :: '[' number ']' | ''
        element   :: symbol isotope count
        separator :: '+' | ' '
        group     :: count element+ | '(' formula ')' count
        grammar   :: group separator formula | group | ''
