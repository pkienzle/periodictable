# This program is public domain

"""
Extensible periodic table of elements

The periodictable package contains mass for the isotopes and density for the
elements.  It calculates xray and neutron scattering information for
isotopes and elements.  Composite values can be calculated from
chemical formula and density.

The table is extensible.  See the user manual for details.

----

Disclaimer:

This data has been compiled from a variety of sources for the user's
convenience and does not represent a critical evaluation by the authors.
While we have made efforts to verify that the values we use match
published values, the values themselves are based on measurements
whose conditions may differ from those of your experiment.

----

"""
__docformat__ = 'restructuredtext en'
__all__ = ['elements', 'neutron_sld','xray_sld','formula'] # and all elements

from . import core
from . import mass
from . import density

# Make a common copy of the table for everyone to use --- equivalent to
# a singleton without incurring any complexity.
elements = core.PeriodicTable()

# Always make mass and density available
mass.init(elements)
density.init(elements)

# Export variables for each element name and symbol.
__all__ += core.define_elements(elements, globals())


def _load_covalent_radius():
    """
    Add covalent_radius property to the elements.

    Note: covalent radii data source is unknown.
    """
    import covalent_radius
    covalent_radius.init(elements)
core.delayed_load(['covalent_radius','covalent_radius_units',
                   'covalent_radius_uncertainty'],
                  _load_covalent_radius)

def _load_crystal_structure():
    """
    Add crystal_structure property to the elements.

    Ashcroft and Mermin
    """
    import crystal_structure
    crystal_structure.init(elements)
core.delayed_load(['crystal_structure'], _load_crystal_structure)

# Delayed loading of neutron properties
def _load_neutron():
    """
    Neutron scattering factors, nuclear_spin and abundance
    properties for elements and isotopes.

    H. Rauch and W. Waschkowski
    ILL Neutron Data Booklet.
    """
    import nsf
    nsf.init(elements)
core.delayed_load(['neutron'],_load_neutron)

# Delayed loading of xray properties
def _load_xray():
    """
    X-ray scattering properties for the elements.

    Center for X-Ray Optics.
    L. Henke, E. M. Gullikson, and J. C. Davis

    K_alpha and K_beta1 emission lines for selected elements.

    D. C. Creagh and J. H. Hubbell
    International Tables for Crystallography Volume C.
    """
    import xsf
    xsf.init(elements)
    xsf.init_spectral_lines(elements)
core.delayed_load(['xray','K_alpha','K_beta1','K_alpha_units','K_beta1_units'],
                  _load_xray)


# Constructors and functions
def formula(value=None, density=None, name=None):
    """
    Chemical formula representation.

    Example initializers:

       string:
          m = formula( "CaCO3+6H2O" )
       sequence of fragments:
          m = formula( [(1,Ca),(2,C),(3,O),(6,[(2,H),(1,O)]] )
       molecular math:
          m = formula( "CaCO3" ) + 6*formula( "H2O" )
       another formula (makes a copy):
          m = formula( formula("CaCO3+6H2O") )
       an atom:
          m = formula( Ca )
       nothing:
          m = formula()

    Additional information can be provided:

       density (g / cm**3)   material density
       name (string) common name for the molecule

    Operations:
       m.atoms returns a dictionary of isotope: count for the
          entire molecule

    Formula strings consist of counts and atoms such as "CaCO3+6H2O".
    Groups can be separated by '+' or space, so "CaCO3 6H2O" works as well.
    Groups and be defined using parentheses, such as "CaCO3(H2O)6".
    Parentheses can nest: "(CaCO3(H2O)6)1"
    Isotopes are represented by index, e.g., "CaCO[18]3+6H2O".
    Counts can be integer or decimal, e.g. "CaCO3+(3HO0.5)2".

    For full details see help(periodictable.formulas.formula_grammar)

    This is designed for calculating molar mass and scattering
    length density, not for representing bonds or atom positions.
    We do preserve the structure of the formula so that it can
    be used as a basis for a rich text representation such as
    matplotlib TeX markup.
    """
    import formulas
    return formulas.Formula(value=value,density=density, name=name)

def neutron_sld(formula,density=None,wavelength=1):
    """
    Compute neutron scattering length densities for molecules.
    Returns the scattering length density, the absorption and
    the incoherent scattering in units of 10**-6 Nb.
    """
    import nsf
    return nsf.neutron_sld(formula,density,wavelength)

def xray_sld(formula,density=None,wavelength=None,energy=None):
    """
    Compute neutron scattering length densities for molecules.
    Returns the scattering length density, the absorption and
    the incoherent scattering in units of 10**-6 Nb.

    Either supply the wavelength (A) or the energy (keV) of the X-rays.
    """
    import xsf
    return xsf.xray_sld(formula,density=density,
                        wavelength=wavelength,energy=energy)

#del core, mass, density
