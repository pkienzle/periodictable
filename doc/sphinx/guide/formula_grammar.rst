.. _formula:


********************
Chemical Composition
********************

Some properties are available for groups of elements.  Groups are specified
as a chemical formula string and either density or cell volume for the crystal
structure.  While it does not provide any information about molecular 
structure, a formula does all com\ *pl*\ ete control over chemical composition. 

* Individual atoms are represented by periodic table symbol.  These are
  case sensitive, so "CO" is different from "Co".
  
* Formula strings consist of counts and atoms such as "CaCO3+6H2O".

* Groups can be separated by '+' or space, so "CaCO3 + 6H2O" works as well. 

* Groups can be defined using parentheses, such as "CaCO3 (H2O)6". 

* Parentheses can nest, e.g., in polyethylene glycol: "HO ((CH2)2O)6 H".

* Isotopes are represented by index, e.g., "CaCO[18]3+6H2O". 

* Counts can be integer or decimal, e.g. "CaCO3+(3HO1.5)2".

A formula string is translated into a formula using 
:mod:`formula class <periodictable.formulas>`. Once the formula has been formed,
you can perform algebra on the entire formula, such as adding
together two formulas to make a more complex compound.

The following is an example of hydrated quartz:

    >>> import periodictable
    >>> SiO2 = periodictable.formula('SiO2')
    >>> hydrated = SiO2 + periodictable.formula('3H2O')
    >>> print hydrated,'mass',hydrated.mass
    SiO2(H2O)3 mass 114.13014
    >>> rho,mu,inc = periodictable.neutron_sld('SiO2+3H2O',density=1.5,wavelength=4.75)
    >>> print hydrated,'neutron sld','%.3g'%rho
    SiO2(H2O)3 neutron sld 0.849
    >>> rho,mu = periodictable.xray_sld(hydrated,density=1.5,
    ... wavelength=periodictable.Cu.K_alpha)
    >>> print hydrated,'X-ray sld','%.3g'%rho
    SiO2(H2O)3 X-ray sld 13.5

Density can be specified directly when the formula is created, or updated
within a formula.  For isotope specific formulas, the density can be given
either as the density of the formula using naturally occurring abundance
if the unit cell is approximately the same, or using the density specific
to those isotopes used.

This makes heavy water density easily specified as:

    >>> D2O = periodictable.formula('D2O',natural_density=1)
    >>> print D2O,"%.4g"%D2O.density
    D2O 1.112

Mixtures can be created by weight or volume ratios, with the density of
the result computed from the density of the materials.  For example, the
following is a 2:1 mixture of water and heavy water:

    >>> H2O = periodictable.formula('H2O',natural_density=1)
    >>> D2O = periodictable.formula('D2O',natural_density=1)
    >>> mix = periodictable.mix_by_volume(H2O,2,D2O,1)
    >>> print mix,"%.4g"%mix.density
    (H2O)2D2O 1.037

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
