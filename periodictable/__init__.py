# This program is public domain
# Author: Paul Kienzle

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
__all__ = ['elements', 'neutron_sld', 'xray_sld',
           'formula', 'mix_by_weight', 'mix_by_volume'] # and all elements
__version__ = "1.6.1"

from . import core
from . import mass
from . import density

# Always make mass and density available
elements = core.PUBLIC_TABLE
mass.init(elements)
density.init(elements)


# Data needed for setup.py when bundling the package into an exe
def data_files():
    """
    Return the data files associated with all periodic table attributes.

    The format is a list of (directory, [files...]) pairs which can be
    used directly in setup(..., data_files=...) for setup.py.

    """
    import os
    import glob
    def _finddata(ext, patterns):
        files = []
        path = core.get_data_path(ext)
        for p in patterns:
            files += glob.glob(os.path.join(path, p))
        return files

    files = [('periodictable-data/xsf',
              _finddata('xsf', ['*.nff', 'read.me', 'f0_WaasKirf.dat'])),
             ('periodictable-data', _finddata('.', ['activation.dat']))]
    return files

# Export variables for each element name and symbol.
__all__ += core.define_elements(elements, globals())

def _load_covalent_radius():
    """
    covalent radius: average atomic radius when bonded to C, N or O.
    """
    from . import covalent_radius
    covalent_radius.init(elements)
core.delayed_load(['covalent_radius',
                   'covalent_radius_units',
                   'covalent_radius_uncertainty'],
                  _load_covalent_radius)

def _load_crystal_structure():
    """
    Add crystal_structure property to the elements.

    Reference:
        *Ashcroft and Mermin.*
    """

    from . import crystal_structure
    crystal_structure.init(elements)
core.delayed_load(['crystal_structure'], _load_crystal_structure)

def _load_neutron():
    """
    Neutron scattering factors, *nuclear_spin* and *abundance*
    properties for elements and isotopes.

    Reference:
        *Rauch. H. and Waschkowski. W., ILL Nuetron Data Booklet.*
    """

    from . import nsf
    nsf.init(elements)
core.delayed_load(['neutron'], _load_neutron, isotope=True)

def _load_neutron_activation():
    """
    Neutron activation calculations for isotopes and formulas.

    Reference:
        *IAEA 273: Handbook on Nuclear Activation Data.*
        *NBSIR 85-3151: Compendium of Benchmark Neutron Field.*
    """
    from . import activation
    activation.init(elements)
core.delayed_load(['neutron_activation'], _load_neutron_activation,
                  element=False, isotope=True)

def _load_xray():
    """
    X-ray scattering properties for the elements.

    Reference:
        *Center for X-Ray optics. Henke. L., Gullikson. E. M., and Davis. J. C.*
    """

    from . import xsf
    xsf.init(elements)
core.delayed_load(['xray'], _load_xray, ion=True)

def _load_emission_lines():
    """
    X-ray emission lines for various elements, including Ag, Pd, Rh, Mo,
    Zn, Cu, Ni, Co, Fe, Mn, Cr and Ti. *K_alpha* is the average of
    K_alpha1 and K_alpha2 lines.
    """

    from . import xsf
    xsf.init_spectral_lines(elements)
core.delayed_load(['K_alpha', 'K_beta1', 'K_alpha_units', 'K_beta1_units'],
                  _load_emission_lines)

def _load_magnetic_ff():
    """
    Magnetic Form Fators. These values are directly from CrysFML.

    Reference:
        *Brown. P. J.(Section 4.4.5)
        International Tables for Crystallography Volume C, Wilson. A.J.C.(ed).*
    """

    from . import magnetic_ff
    magnetic_ff.init(elements)
core.delayed_load(['magnetic_ff'], _load_magnetic_ff)


# Constructors and functions
def formula(*args, **kw):
    """
    Chemical formula representation.

    Example initializers:

       string:
          m = formula( "CaCO3+6H2O" )
       sequence of fragments:
          m = formula( [(1, Ca), (2, C), (3, O), (6, [(2, H), (1, O)]] )
       molecular math:
          m = formula( "CaCO3" ) + 6*formula( "H2O" )
       another formula (makes a copy):
          m = formula( formula("CaCO3+6H2O") )
       an atom:
          m = formula( Ca )
       nothing:
          m = formula()

    Additional information can be provided:

       density (|g/cm^3|)   material density
       natural_density (|g/cm^3|) material density with natural abundance
       name (string) common name for the molecule
       table (PeriodicTable) periodic table with customized data

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

    The chemical formula is designed for simple calculations such
    as molar mass, not for representing bonds or atom positions.
    However, we preserve the structure of the formula so that it can
    be used as a basis for a rich text representation such as
    matplotlib TeX markup.
    """
    from . import formulas
    return formulas.formula(*args, **kw)

def mix_by_weight(*args, **kw):
    """
    Generate a mixture which apportions each formula by weight.

    :Parameters:

        *formula1* : Formula OR string
            Material

        *quantity1* : float
            Relative quantity of that material

        *formula2* : Formula OR string
            Material

        *quantity2* : float
            Relative quantity of that material

        ...

        *density* : float
            Density of the mixture, if known

        *natural_density* : float
            Density of the mixture with natural abundances, if known.

        *name* : string
            Name of the mixture

    :Returns:

        *formula* : Formula

    If density is not given, then it will be computed from the density
    of the components, assuming equal volume.
    """
    from . import formulas
    return formulas.mix_by_weight(*args, **kw)

def mix_by_volume(*args, **kw):
    """
    Generate a mixture which apportions each formula by volume.

    :Parameters:

        *formula1* : Formula OR string
            Material

        *quantity1* : float
            Relative quantity of that material

        *formula2* : Formula OR string
            Material

        *quantity2* : float
            Relative quantity of that material

        ...

        *density* : float
            Density of the mixture, if known

        *natural_density* : float
            Density of the mixture with natural abundances, if known.

        *name* : string
            Name of the mixture

    :Returns:

        *formula* : Formula

    Densities are required for each of the components.  If the density of
    the result is not given, it will be computed from the components
    assuming the components take up no more nor less space because they
    are in the mixture.
    """
    from . import formulas
    return formulas.mix_by_volume(*args, **kw)


def neutron_sld(*args, **kw):
    """
    Compute neutron scattering length densities for molecules.

    Returns scattering length density (real, imaginary and incoherent).

    See :class:`periodictable.nsf.neutron_sld` for details.
    """
    from . import nsf
    return nsf.neutron_sld(*args, **kw)

def neutron_scattering(*args, **kw):
    """
    Compute neutron scattering cross sections for molecules.

    Returns scattering length density (real, imaginary and incoherent),
    cross sections (coherent, absorption, incoherent) and penetration
    depth.

    See :func:`periodictable.nsf.neutron_scattering` for details.
    """
    from . import nsf
    return nsf.neutron_scattering(*args, **kw)

def xray_sld(*args, **kw):
    """
    Compute neutron scattering length densities for molecules.

    Either supply the wavelength (A) or the energy (keV) of the X-rays.

    Returns scattering length density (real, imaginary).

    See :class:`periodictable.xsf.Xray` for details.
    """
    from . import xsf
    return xsf.xray_sld(*args, **kw)


#del core, mass, density
