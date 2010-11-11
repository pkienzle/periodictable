.. _formula:


********************
Chemical Composition
********************

Some properties are available for groups of elements.  Groups are specified
as a chemical formula string and either density or cell volume for the crystal
structure.  While it does not provide any information about molecular 
structure, a formula does provide complete control over chemical composition. 

A formula string is translated into a formula using 
:func:`periodictable.formulas.formula`:

* Formula strings consist of counts and atoms, where individual atoms are 
  represented by periodic table symbol.  The atoms are case sensitive, 
  so "CO" is different from "Co".  Here is an example of calcium carbonate:
 
    >>> from periodictable import formula
    >>> print formula("CaCO3")
    CaCO3
  
* Formulas can contain multiple groups separated by space or plus or by using
  parentheses.  Whole groups can have a repeat count.  The following are
  equivalent definitions of hydrated calcium carbonate:
  
    >>> print formula("CaCO3+6H2O")
    CaCO3(H2O)6
    >>> print formula("CaCO3 6H2O")
    CaCO3(H2O)6
    >>> print formula("CaCO3(H2O)6")
    CaCO3(H2O)6

* Parentheses can nest, e.g., in polyethylene glycol:

    >>> print formula("HO ((CH2)2O)6 H")
    HO((CH2)2O)6H

* Isotopes are represented by index, such as O[18] = :sup:`18`\ O:

    >>> print formula("CaCO[18]3+6H2O")
    CaCO[18]3(H2O)6

* Counts can be integer or decimal:

    >>> print formula("CaCO3+(3HO1.5)2")
    CaCO3((HO1.5)3)2

* Empty formulas are supported, e.g., for air or vacuum:
    
    >>> print formula()
    <BLANKLINE>
    >>> formula()
    formula('')

The grammar used for parsing formula strings is the following:

::

    number    :: [1-9][0-9]*
    fraction  :: (number | [0] | nothing) '.' [0-9]*
    count     :: number | fraction | nothing
    symbol    :: [A-Z][a-z]*
    isotope   :: '[' number ']' | nothing
    element   :: symbol isotope count
    separator :: '+' | nothing
    group     :: count element+ | '(' formula ')' count
    grammar   :: group space separator space formula | group | nothing


Formulas can also be constructed from atoms or other formulas:

* A simple formula can be created from a bare atom:

    >>> from periodictable import Ca, C, O, H
    >>> print formula(Ca)
    Ca

* More complex structures will require a sequences of counts and fragments.
  The fragment itself can be a structure:

    >>> print formula( [ (1,Ca), (1,C), (3,O), (6,[(2,H),(1,O)]) ] )
    CaCO3(H2O)6

* Structures can also be built with simple formula math:
    
    >>> print formula("CaCO3") + 6*formula("H2O")
    CaCO3(H2O)6

* Formulas can be easily cloned:
    
    >>> print formula( formula("CaCO3+6H2O"))
    CaCO3(H2O)6

Density
-------

Density can be specified directly when the formula is created, or updated
within a formula.  For isotope specific formulas, the density can be given
either as the density of the formula using naturally occurring abundance
if the unit cell is approximately the same, or using the density specific
to those isotopes used.

This makes heavy water density easily specified as:

    >>> D2O = formula('D2O',natural_density=1)
    >>> print D2O,"%.4g"%D2O.density
    D2O 1.112

Density can also be estimated from the volume of the unit cell, either
by using the covalent radii of the constituent atoms and assuming some
packing factor, or by knowing the lattice parameters of the crystal
which makes up the material.  Standard packing factors for hcp, fcc,
bcc, cubic and diamond on uniform spheres can be used if the components
are of about the same size.  The formula should specify the number of
atoms in the unit cell, which is 1 for cubic, 2 for bcc and 4 for fcc.  
Be sure to use the molecular mass (M.molecular_mass in g) rather 
than the molar mass (M.mass in u = g/mol) in your calculations.

Because the packing fraction method relies on the covalent radius
estimate it is not very accurate:

    >>> from periodictable import elements, formula
    >>> Fe = formula("2Fe")  # bcc lattice has 2 atoms per unit cell
    >>> Fe.density = Fe.molecular_mass/Fe.volume('bcc')
    >>> print "%.3g"%Fe.density
    6.55
    >>> print "%.3g"%elements.Fe.density
    7.87

Using lattice parameters the results are much better:

    >>> Fe.density = Fe.molecular_mass/Fe.volume(a=2.8664)
    >>> print "%.3g"%Fe.density
    7.88
    >>> print "%.3g"%elements.Fe.density
    7.87

Mixtures
--------

Mixtures can be created by weight or volume ratios, with the density of
the result computed from the density of the materials.  For example, the
following is a 2:1 mixture of water and heavy water:

    >>> from periodictable import formula, mix_by_volume, mix_by_weight
    >>> H2O = formula('H2O',natural_density=1)
    >>> D2O = formula('D2O',natural_density=1)
    >>> mix = mix_by_volume(H2O,2,D2O,1)
    >>> print mix,"%.4g"%mix.density
    (H2O)2D2O 1.037
    
Note that this is different from a 2:1 mixture by weight:

    >>> mix = mix_by_weight(H2O,2,D2O,1)
    >>> print mix,"%.4g"%mix.density
    (H2O)2.2234D2O 1.035

Derived values
--------------

Once a formula has been created, it can be used for summary calculations.
The following is an example of hydrated quartz, which shows how to
compute molar mass and neutron/xray scattering length density:

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
