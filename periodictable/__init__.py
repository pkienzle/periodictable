# This program is public domain

"""
Extensible periodic table of elements

The periodictable package contains mass for the isotopes and density for the
elements. It calculates xray and neutron scattering information for
isotopes and elements. Composite values can be calculated from
chemical formula and density.

The table is extensible. See the user manual for details.

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

# Data needed for setup.py when bundling the package into an exe
def data_files():
    """
    Return the data files associated with all periodic table attributes.
    
    The format is a list of (directory, [files...]) pairs which can be
    used directly in setup(...,data_files=...) for setup.py.

    """
    import os, glob
    def _finddata(ext, patterns):
        files = []
        path = core.get_data_path(ext)
        for p in patterns:
            files += glob.glob(os.path.join(path,p))
            return files

    data_files = [('periodictable-data/xsf', 
                   _finddata('xsf', ['*.nff','read.me','f0_WaasKirf.dat']))]
    return data_files

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

    .. Note:: covalent radii data source is unknown.
    """
    import covalent_radius
    covalent_radius.init(elements)
core.delayed_load(['covalent_radius'],
                  _load_covalent_radius)

def _load_covalent_radius():
    """
    Add covalent_radius units property to the elements.
    """

    import covalent_radius
    covalent_radius.init(elements)
core.delayed_load(['covalent_radius_units'],
                  _load_covalent_radius)

def _load_covalent_radius():
    """
    Add covalent_radius uncertainty property to the elements.
    """

    import covalent_radius
    covalent_radius.init(elements)
core.delayed_load(['covalent_radius_uncertainty'],
                  _load_covalent_radius)

def _load_crystal_structure():
    """
    Add crystal_structure property to the elements.

    Reference:
        *Ashcroft and Mermin.*
    """

    import crystal_structure
    crystal_structure.init(elements)
core.delayed_load(['crystal_structure'], _load_crystal_structure)

def _load_neutron():
    """
    Neutron scattering factors, *nuclear_spin* and *abundance*
    properties for elements and isotopes.

    Reference:
        *Rauch. H. and Waschkowski. W., ILL Nuetron Data Booklet.*
    """

    import nsf
    nsf.init(elements)
core.delayed_load(['neutron'],_load_neutron, isotope=True)

def _load_xray():
    """
    X-ray scattering properties for the elements.
     
    Reference:
        *Center for X-Ray optics. Henke. L., Gullikson. E. M., and Davis. J. C.*
    """

    import xsf
    xsf.init(elements)
core.delayed_load(['xray'], _load_xray, ion=True)

def _load_emission_lines():
    """
    X-ray emission lines for various elements, including Ag, Pd, Rh, Mo,
    Zn, Cu, Ni, Co, Fe, Mn, Cr and Ti. *K_alpha* is the average of
    K_alpha1 and K_alpha2 lines.
    """
 
    import xsf
    xsf.init_spectral_lines(elements)
core.delayed_load(['K_alpha','K_beta1','K_alpha_units','K_beta1_units'],
                  _load_emission_lines)

def _load_magnetic_ff():
    """
    Magnetic Form Fators. These values are directly from CrysFML.

    Reference:
        *Brown. P. J.(Section 4.4.5) International Tables for Crystallography Volume C, Wilson. A.J.C.(ed).*
    """

    import magnetic_ff
    magnetic_ff.init(elements)
core.delayed_load(['magnetic_ff'], _load_magnetic_ff)

def _load_ionic_radius():
    """
    Ionic radii for various charges.These values are directly from CrysFML. 
    """
    import ionic_radius
    ionic_radius.init(elements)
core.delayed_load(['ionic_radius'], _load_ionic_radius)



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

    Returns scattering length density (real, imaginary and incoherent).
    
    See :class:`periodictable.nsf.Neutron` for details.  
    """
    import nsf
    return nsf.neutron_sld(formula,density,wavelength)

def xray_sld(formula,density=None,wavelength=None,energy=None):
    """
    Compute neutron scattering length densities for molecules.
    
    Either supply the wavelength (A) or the energy (keV) of the X-rays.

    Returns scattering length density (real, imaginary).
    
    See :class:`periodictable.xsf.Xray` for details.
    """
    import xsf
    return xsf.xray_sld(formula,density=density,
                        wavelength=wavelength,energy=energy)


#del core, mass, density
