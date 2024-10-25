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
__version__ = "2.0.0-pre"

__all__ = ['elements'] # Lazy symbols and individual elements added later

import importlib

from . import core
from . import mass
from . import density

_LAZY_MODULES = [
    'activation',
    'covalent_radius',
    'crystal_structure',
    'nsf',
    'xsf',
    'magnetic_ff',
    ]
_LAZY_LOAD = {
    'formula': 'formulas',
    'mix_by_weight': 'formulas',
    'mix_by_volume': 'formulas',
    'neutron_sld': 'nsf',
    'neutron_scattering': 'nsf',
    'xray_sld': 'xsf',
}
def __getattr__(name: str):
    """
    Lazy loading of modules and symbols from other modules. This is
    equivalent to using "from .formulas import formula" etc in __init__
    except that the import doesn't happen until the symbol is referenced.
    Using "from periodictable import formula" will import the symbol immediately.
    "from periodictable import *" will import all symbols, including the lazy
    """
    module_name = _LAZY_LOAD.get(name, None)
    if module_name is not None:
        # Lazy symbol: fetch name from the target module
        #print(f"from {__name__}.{module_name} import {name} [lazy]")
        module = importlib.import_module(f'{__name__}.{module_name}')
        symbol = getattr(module, name)
        globals()[name] = symbol
        return symbol
    if name in _LAZY_MODULES:
        # Lazy module: just need to import it
        #print(f"import {__name__}.{name} [lazy]")
        return importlib.import_module(f'{__name__}.{name}')
    raise AttributeError(f"module '{__name__}' has not attribute '{name}'")
def __dir__():
    return __all__
# Support 'from periodictable import *' and 'dir(periodictable)'
__all__ = [*__all__, *_LAZY_MODULES, *_LAZY_LOAD.keys()]

# Always make mass and density available
elements = core.PUBLIC_TABLE
mass.init(elements)
density.init(elements)
del mass, density

# Add element name and symbol (e.g. nickel and Ni) to the public attributes.
__all__ += core.define_elements(elements, globals())

# Lazy loading of element and isotope attributes, e.g., Ni.covalent_radius
core.delayed_load(
    ['covalent_radius', 'covalent_radius_units', 'covalent_radius_uncertainty'],
    lambda: covalent_radius.init(elements))
core.delayed_load(
    ['crystal_structure'],
    lambda: crystal_structure.init(elements))
core.delayed_load(
    ['neutron'],
    lambda: nsf.init(elements), isotope=True)
core.delayed_load(
    ['neutron_activation'],
    lambda: activation.init(elements),
    element=False, isotope=True)
core.delayed_load(
    ['xray'],
    lambda: xsf.init(elements),
    ion=True)
core.delayed_load(
    ['K_alpha', 'K_beta1', 'K_alpha_units', 'K_beta1_units'],
    lambda: xsf.init_spectral_lines(elements))
core.delayed_load(
    ['magnetic_ff'],
    lambda: magnetic_ff.init(elements))

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
              _finddata('xsf', ['*.nff', 'read.me'])),
             ('periodictable-data', _finddata('.', ['activation.dat', 'f0_WaasKirf.dat']))]
    return files
