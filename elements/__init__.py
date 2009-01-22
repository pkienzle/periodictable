# This program is public domain

"""
Information about elements and isotopes.

The elements package contains mass for the isotopes and density for the
elements.  It calculates xray and neutron scattering information for
isotopes, elements and molecules.

The table is extensible.  See help(elements.elements) for details.

Basic usage
===========

Access particular elements by name::

    from elements import hydrogen
    print "H mass", hydrogen.mass, hydrogen.mass_units

Access particular elements as symbols::

    from elements import H,B,Cu,Ni
    print "H mass", H.mass, H.mass_units
    print "B absorption", B.neutron.absorption, B.neutron.absorption_units
    print "Ni f1/f2 for Cu K-alpha X-rays", Ni.xray.scattering_factors(Cu.K_alpha)

Access isotopes using mass number subscripts::

    print "58-Ni vs 62-Ni scattering", \\
        Ni[58].neutron.coherent, Ni[62].neutron.coherent

Access elements indirectly::

    import elements
    print "Cd density", elements.Cd.density, elements.Cd.density_units

Import all elements::

    from elements import *

Deuterium and tritium are special isotopes named D and T
some neutron information is available as 'n'::

    print "D mass",D.mass
    print "neutron mass",n.mass

Process all the elements::

    for el in elements.table: # if you used "import elements"
        print el.symbol,el.name

    for el in periodic_table: # if you used "from elements import *"
        print el.symbol,el.number

Process all the isotopes for an element::

    for iso in elements.Fe:
        print iso,iso.mass

Missing properties generally evaluate to None::

    print "Radon density",elements.Rn.density

Helper function for listing only the defined properties::

    periodic_table.list('symbol','K_alpha',format="%s K-alpha=%s")

Work with molecules::

    SiO2 = elements.molecule('SiO2')
    hydrated = SiO2 + elements.molecule('3H2O')
    print hydrated,'mass',hydrated.mass
    rho,mu,inc = elements.neutron_sld('SiO2+3H2O',density=1.5,wavelength=4.75)
    print hydrated,'neutron sld','%.3g'%rho
    rho,mu = elements.xray_sld(hydrated,density=1.5,wavelength=Cu.K_alpha)
    print hydrated,'X-ray sld','%.3g'%rho

-----------
Disclaimer:

This data has been compiled from the above sources for the user's 
convenience and does not represent a critical evaluation by the authors.
While we have made efforts to verify that the values we use match
published values, the values themselves are based on measurements
whose conditions may differ from those of your experiment.
----
"""
__docformat__ = 'restructuredtext en'

# Pull in periodic table and elements
import core
import mass
import density
from core import *
__all__ = core.__all__ + ['neutron_sld','xray_sld','molecule']

# Allow elements.table as a shorthand for elements.periodic_table
table = periodic_table

def _load_covalent_radius():
    """
    Add covalent_radius property to the elements.
    
    Note: covalent radii data source is unknown.
    """
    import covalent_radius
    covalent_radius._init()
core.delayed_load(['covalent_radius'],_load_covalent_radius)

def _load_crystal_structure():
    """
    Add crystal_structure property to the elements.
    
    Ashcroft and Mermin
    """
    import crystal_structure
core.delayed_load(['crystal_structure'],_load_covalent_radius)

# Delayed loading of neutron properties
def _load_neutron():
    """
    Neutron scattering factors, nuclear_spin and abundance
    properties for elements and isotopes.

    H. Rauch and W. Waschkowski
    ILL Neutron Data Booklet.
    """
    import nsf
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
core.delayed_load(['xray','K_alpha','K_beta1'],_load_xray)
core.Element.K_alpha_units = "angstrom"
core.Element.K_beta1_units = "angstrom"


# Constructors and functions
def molecule(value=None, density=None, name=None):
    """
    Molecule representation.  

    Example initializers:

       string: 
          m = Molecule( "CaCO3+6H2O" )
       sequence of fragments:
          m = Molecule( [(1,Ca),(2,C),(3,O),(6,[(2,H),(1,O)]] )
       molecular math:
          m = Molecule( "CaCO3" ) + 6*Molecule( "H2O" )
       another molecule (makes a copy):
          m = Molecule( Molecule("CaCO3+6H2O") )
       an atom:
          m = Molecule( Ca )
       nothing:
          m = Molecule()
          
    Additional information can be provided:
    
       density (g / cm**3)   material density
       name (string) common name for the molecule

    Operations:
       m.atoms returns a dictionary of isotope: count for the
          entire molecule

    Molecule strings consist of counts and atoms such as "CaCO3+6H2O".  
    Groups can be separated by '+' or space, so "CaCO3 6H2O" works as well. 
    Groups and be defined using parentheses, such as "CaCO3(H2O)6".
    Parentheses can nest: "(CaCO3(H2O)6)1"
    Isotopes are represented by index, e.g., "CaCO[18]3+6H2O". 
    Counts can be integer or decimal, e.g. "CaCO3+(3HO0.5)2".

    For full details see help(elements.molecules.molecule_grammar)
        
    This is designed for calculating molar mass and scattering 
    length density, not for representing bonds or atom positions.  
    We do preserve the structure of the formula so that it can 
    be used as a basis for a rich text representation such as 
    matplotlib TeX markup.
    """
    import molecules
    return molecules.Molecule(value=value,density=density, name=name)

def neutron_sld(molecule,density=None,wavelength=1):
    """
    Compute neutron scattering length densities for molecules.
    Returns the scattering length density, the absorption and
    the incoherent scattering in units of 10**-6 Nb.
    """
    import nsf
    return nsf.neutron_sld(molecule,density,wavelength)

def xray_sld(molecule,density=None,wavelength=None,energy=None):
    """
    Compute neutron scattering length densities for molecules.
    Returns the scattering length density, the absorption and
    the incoherent scattering in units of 10**-6 Nb.

    Either supply the wavelength (A) or the energy (keV) of the X-rays.
    """
    import xsf
    return xsf.xray_sld(molecule,density=density,
                        wavelength=wavelength,energy=energy)

del core, mass, density

