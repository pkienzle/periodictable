# This program is public domain
"""
Periodic table definition, with classes for elements and isotopes.

PeriodTable
   Creates an empty periodic table.  Properties are added to it
   separately.  Note this is not a singleton.  But unless you are
   very careful you should use the default table defined for the
   package.

Element
   A class to hold the properties for each element.

Isotope
   A class to hold the properties for individual isotopes.


Extending the table
-------------------

The periodic table is extensible.  Third party packages can
add attributes to the table, and they will appear in all of
the elements.  For example:

    import danse.gammaray  # loads gammaray data
    from periodictable import *  # loads symbols for H, He, ...

    for el in periodictable.elements: print el.symbol, el.gammadata

To implement this, you will need the following in gammaray/core.py:

    from periodictable.core import elements, Element

    def _init():
        if 'gammadata' in elements.properties: return
        elements.properties.append('gammadata')

        # Set the default, if any
        Element.gamma = None

        # Set the units
        Element.gamma_units = "Ev"

        # Load the data
        for s,data in gamma_table.iteritems():
            el = elements.symbol(s)
            el.gammadata = data

    # Define the data
    gamma_table = dict(
        Si="Silicon gamma values",
        O="Oxygen gamma values",
        )

    _init()

You can use similar tricks for isotope specific data:

    from periodictable.core import elements, Isotope

    def _init():
        if 'shells' in elements.properties: return
        elements.properties.append('shells')

        # Set the default.  This is required, even if it is only
        # setting it to None.  If the attribute is missing then the
        # isotope data reverts to the element to supply the value,
        # which is almost certainly not what you want.
        Isotope.shells = None

        # Load the data
        for symbol,data in shell_table.iteritems():
            el = elements.symbol(symbol)
            for iso,isodata in data.iteritems():
                el[iso].shells = isodata

    # Define the data
    shell_table = dict(
        Fe={56: "56-Fe shell info",
            58: "58-Fe shell info",
            }
        )

    _init()

Since your data table is in its own package, you need an __init__.py to
create the initial periodic table.  It should contain something like:

     import periodictable
     ...
     del periodictable  # Clean up symbol space

You can also define attributes on import that are not loaded directly.
For example, if you don't want to load all the isotope information for
shells immediately, then you can use delayed_load in __init__.py:


     from periodictable.core import delayed_load

     # Delayed loading of shell info
     def _load_shell():
         '''
         Electron shell information for isotopes.

         T. Student, Tables of Shell Information
         '''
         import shelltable
     delayed_load(['shells'],_load_shell)

The first argument to delayed_load is the list of all attributes that will
be defined when the module is loaded.  The second argument is the loading
function, whose docstring will appear as the attribute description for
each attribute in the first list.
"""
import copy

# Note: __all__ will include all the elements and elements; it is
# defined below.

def delayed_load(all_props,loader,element=True,isotope=False):
    """
    Delayed loading of an element property table.  When any of props
    is first accessed the loader will be called to load the associated
    data.  The help string starts out as the help string for the loader
    function.  If it is an isotope property, be sure to set the
    keyword isotope=True.
    """
    def clearprops():
        """
        Remove the properties so that the attribute can be accessed
        directly.
        """
        if element:
            for p in all_props:
                delattr(Element, p)
        if isotope:
            for p in all_props:
                delattr(Isotope, p)

    def getter(propname):
        """
        Property getter for attribute propname.

        The first time the prop is accessed, the prop itself will be
        deleted and the data loader for the property will be called
        to set the real values.  Subsequent references to the property
        will be to the actual data.
        """
        def getfn(el):
            #print "get",el,propname
            clearprops()
            loader()
            return getattr(el, propname)
        return getfn

    def setter(propname):
        """
        Property setter for attribute propname.

        This function is assumed to be called when the data loader for the
        attribute is called before the property is referenced (for example,
        if somebody imports periodictable.xsf before referencing Ni.xray). In
        this case, we simply need to clear the delayed load property and
        let the loader set the values as usual.

        If the user tries to override a value in the table before first
        referencing the table, then the above assumption is false.  E.g.,
        "Ni.K_alpha=5" followed by "print Cu.K_alpha" will yield an
        undefined Cu.K_alpha.  This will be difficult for future users
        to debug.
        """
        def setfn(el, value):
            #print "set",el,propname,value
            clearprops()
            setattr(el, propname, value)
        return setfn

    if element:
        for p in all_props:
            prop = property(getter(p), setter(p), doc=loader.__doc__)
            setattr(Element, p, prop)

    if isotope:
        for p in all_props:
            prop = property(getter(p), setter(p), doc=loader.__doc__)
            setattr(Isotope, p, prop)


# Define the element names from the element table.
class _PeriodicTable(object):
    """
    Defines the period table of the elements with isotopes.

    Individidual elements are accessed by name, symbol or atomic number.
    Individual isotopes are addressable by element[mass_number] or
    elements.isotope('#-Xx').

    For example, the following all retrieve iron::

        elements[26]
        elements.Fe
        elements.symbol('Fe')
        elements.name('iron')
        elements.isotope('Fe')

    To get iron-56, use::

        elements[26][56]
        elements.Fe[56]
        elements.isotope('56-Fe')

    Deuterium and tritium are defined as D and T.  Some
    neutron properties are available in elements[0]

    To show all the elements in the table, use the iterator:

        for element in elements:  # lists the element symbols
            print el.symbol,el.number

    Properties can be added to the elements as needed, including mass, nuclear
    and X-ray scattering cross sections.  See the help(periodictable)
    for details.
    """
    def __init__(self):
        """
        Create one element for each entry in the periodic table.
        """
        self.properties = []
        self._element = {}
        for Z,name in element_symbols.iteritems():
            element = Element(name[0].lower(),name[1],Z)
            self._element[element.number] = element
            setattr(self,element.symbol,element)

        # There are two specially named isotopes D and T
        self.D = self.H.add_isotope(2)
        self.D.name = 'deuterium'
        self.D.symbol = 'D'
        self.T = self.H.add_isotope(3)
        self.T.name = 'tritium'
        self.T.symbol = 'T'

    def __getitem__(self, Z):
        """
        Retrieve element Z
        """
        return self._element[Z]

    def __iter__(self):
        """
        Process the elements in Z order
        """
        Zlist = self._element.keys()
        Zlist.sort()
        for Z in Zlist:
            yield self._element[Z]

    def symbol(self, input):
        """
        Lookup the symbol in the periodic table (and make sure it is an element
        or deuterium or tritium (D or T).
        """
        if hasattr(self,input):
            value = getattr(self,input)
            if isinstance(value,(Element,Isotope)):
                return value
        raise ValueError("unknown element "+input)

    def name(self, input):
        """
        Lookup the full name of the element in the period table.
        """
        for el in self:
            if input == el.name: return el
        if input == self.D.name: return self.D
        if input == self.T.name: return self.T
        raise ValueError("unknown element "+input)

    def isotope(self, input):
        """
        Lookup the element or isotope in the periodic table. Elements
        are assumed to be given by the standard element symbols.  Isotopes
        are given by number-symbol, or D and T for 2-H and 3-H.
        """
        # Parse #-Sym or Sym
        # If # is not an integer, set isotope to -1 so that the isotope
        # lookup will fail later.
        parts = input.split('-')
        if len(parts) == 1:
            isotope = 0
            symbol = parts[0]
        elif len(parts) == 2:
            try:
                isotope = int(parts[0])
            except:
                isotope = -1
            symbol = parts[1]
        else:
            symbol = ''
            isotope = -1

        # All elements are attributes of the table
        # Check that the attribute is an Element or an Isotope (for D or T)
        # If it is an element, check that the isotope exists
        if hasattr(self,symbol):
            attr = getattr(self,symbol)
            if isinstance(attr,Element):
                # If no isotope, return the element
                if isotope == 0:
                    return attr
                # If isotope, check that it is valid
                if isotope in attr.isotopes:
                    return attr[isotope]
            elif isinstance(attr,Isotope):
                # D,T must not have an associated isotope; D[4] is meaningless.
                if isotope == 0:
                    return attr

        # If we can't parse the string as an element or isotope, raise an error
        raise ValueError("unknown element "+input)

    def list(self,*props,**kw):
        """
        list('prop1','prop2',...,format='format')
            Print a list of elements with the given set of properties.

        If format is given, then it must contain one %s for each property.
        """
        format = kw.pop('format',None)
        assert len(kw) == 0
        for el in self:
            try:
                L = tuple(getattr(el,p) for p in props)
            except AttributeError:
                # Skip elements which don't define all the attributes
                continue

            if format is None:
                print " ".join(str(p) for p in L)
            else:
                try:
                    print format%L
                except:
                    print "format",format,"args",L
                    raise

class Isotope(object):
    """Periodic table entry for an individual isotope.

    An isotope is associated with an element.  In addition to the element
    properties (symbol, name, atomic number), it has specific isotope
    properties (isotope number, nuclear spin, relative abundance).
    Properties not specific to the isotope (e.g., x-ray scattering factors)
    are retrieved from the associated element.
    """
    def __init__(self,element,isotope_number):
        self.element = element
        self.isotope = isotope_number
    def __getattr__(self,attr):
        return getattr(self.element,attr)
    def __str__(self):
        # Deuterium and Tritium are special
        if 'symbol' in self.__dict__:
            return self.symbol
        else:
            return "%d-%s"%(self.isotope,self.element.symbol)
    def __repr__(self):
        return "%s[%d]"%(self.element.symbol,self.isotope)
    def __getstate__(self):
        """
        Can't pickle isotopes without breaking extensibility.
        Suppress it for now.
        """
        raise TypeError("cannot copy or pickle isotopes")

class Element(object):
    """Periodic table entry for an element.

    An element is a name, symbol and number, plus a set of properties.
    Individual isotopes can be referenced as element[isotope_number].
    """
    def __init__(self,name,symbol,Z):
        self.name = name
        self.symbol = symbol
        self.number = Z
        self._isotopes = {} # The actual isotopes
    def _getisotopes(self):
        L = self._isotopes.keys()
        L.sort()
        return L
    isotopes = property(_getisotopes,doc="List of all isotopes")
    def add_isotope(self,number):
        if number not in self._isotopes:
            self._isotopes[number] = Isotope(self,number)
        return self._isotopes[number]
    def __getitem__(self,number):
        return self._isotopes[number]

    def __iter__(self):
        """
        Process the elements in Z order
        """
        for isotope in self.isotopes:
            yield self._isotopes[isotope]

    # Note: using repr rather than str for the element symbol so
    # that lists of elements print nicely.  Since elements are
    # effectively singletons, the symbol name is the representation
    # of the instance.
    def __repr__(self):
        return self.symbol

    def __getstate__(self):
        """
        Can't pickle elements without breaking extensibility.
        Suppress it for now.
        """
        raise TypeError("cannot copy or pickle elements")

element_symbols = {
    0: ['Neutron','n'],
    1: ['Hydrogen', 'H'],
    2: ['Helium', 'He'],
    3: ['Lithium', 'Li'],
    4: ['Beryllium', 'Be'],
    5: ['Boron', 'B'],
    6: ['Carbon', 'C'],
    7: ['Nitrogen', 'N'],
    8: ['Oxygen', 'O'],
    9: ['Fluorine', 'F'],
    10: ['Neon', 'Ne'],
    11: ['Sodium', 'Na'],
    12: ['Magnesium', 'Mg'],
    13: ['Aluminum', 'Al'],
    14: ['Silicon', 'Si'],
    15: ['Phosphorus', 'P'],
    16: ['Sulfur', 'S'],
    17: ['Chlorine', 'Cl'],
    18: ['Argon', 'Ar'],
    19: ['Potassium', 'K'],
    20: ['Calcium', 'Ca'],
    21: ['Scandium', 'Sc'],
    22: ['Titanium', 'Ti'],
    23: ['Vanadium', 'V'],
    24: ['Chromium', 'Cr'],
    25: ['Manganese', 'Mn'],
    26: ['Iron', 'Fe'],
    27: ['Cobalt', 'Co'],
    28: ['Nickel', 'Ni'],
    29: ['Copper', 'Cu'],
    30: ['Zinc', 'Zn'],
    31: ['Gallium', 'Ga'],
    32: ['Germanium', 'Ge'],
    33: ['Arsenic', 'As'],
    34: ['Selenium', 'Se'],
    35: ['Bromine', 'Br'],
    36: ['Krypton', 'Kr'],
    37: ['Rubidium', 'Rb'],
    38: ['Strontium', 'Sr'],
    39: ['Yttrium', 'Y'],
    40: ['Zirconium', 'Zr'],
    41: ['Niobium', 'Nb'],
    42: ['Molybdenum', 'Mo'],
    43: ['Technetium', 'Tc'],
    44: ['Ruthenium', 'Ru'],
    45: ['Rhodium', 'Rh'],
    46: ['Palladium', 'Pd'],
    47: ['Silver', 'Ag'],
    48: ['Cadmium', 'Cd'],
    49: ['Indium', 'In'],
    50: ['Tin', 'Sn'],
    51: ['Antimony', 'Sb'],
    52: ['Tellurium', 'Te'],
    53: ['Iodine', 'I'],
    54: ['Xenon', 'Xe'],
    55: ['Cesium', 'Cs'],
    56: ['Barium', 'Ba'],
    57: ['Lanthanum', 'La'],
    58: ['Cerium', 'Ce'],
    59: ['Praseodymium', 'Pr'],
    60: ['Neodymium', 'Nd'],
    61: ['Promethium', 'Pm'],
    62: ['Samarium', 'Sm'],
    63: ['Europium', 'Eu'],
    64: ['Gadolinium', 'Gd'],
    65: ['Terbium', 'Tb'],
    66: ['Dysprosium', 'Dy'],
    67: ['Holmium', 'Ho'],
    68: ['Erbium', 'Er'],
    69: ['Thulium', 'Tm'],
    70: ['Ytterbium', 'Yb'],
    71: ['Lutetium', 'Lu'],
    72: ['Hafnium', 'Hf'],
    73: ['Tantalum', 'Ta'],
    74: ['Tungsten', 'W'],
    75: ['Rhenium', 'Re'],
    76: ['Osmium', 'Os'],
    77: ['Iridium', 'Ir'],
    78: ['Platinum', 'Pt'],
    79: ['Gold', 'Au'],
    80: ['Mercury', 'Hg'],
    81: ['Thallium', 'Tl'],
    82: ['Lead', 'Pb'],
    83: ['Bismuth', 'Bi'],
    84: ['Polonium', 'Po'],
    85: ['Astatine', 'At'],
    86: ['Radon', 'Rn'],
    87: ['Francium', 'Fr'],
    88: ['Radium', 'Ra'],
    89: ['Actinium', 'Ac'],
    90: ['Thorium', 'Th'],
    91: ['Protactinium', 'Pa'],
    92: ['Uranium', 'U'],
    93: ['Neptunium', 'Np'],
    94: ['Plutonium', 'Pu'],
    95: ['Americium', 'Am'],
    96: ['Curium', 'Cm'],
    97: ['Berkelium', 'Bk'],
    98: ['Californium', 'Cf'],
    99: ['Einsteinium', 'Es'],
    100: ['Fermium', 'Fm'],
    101: ['Mendelevium', 'Md'],
    102: ['Nobelium', 'No'],
    103: ['Lawrencium', 'Lr'],
    104: ['Rutherfordium', 'Rf'],
    105: ['Dubnium', 'Db'],
    106: ['Seaborgium', 'Sg'],
    107: ['Bohrium', 'Bh'],
    108: ['Hassium', 'Hs'],
    109: ['Meitnerium', 'Mt'],
    110: ['Ununnilium', 'Uun'],
    111: ['Unununium', 'Uuu'],
    112: ['Ununbium', 'Uub'],
    114: ['Ununquadium', 'Uuq'],
    116: ['Ununhexium', 'Uuh'],
}


# Make a common copy of the table for everyone to use --- equivalent to
# a singleton without incurring any complexity; this must be done after
# the table data is defined.
elements = _PeriodicTable()
for el in elements:
    exec el.symbol+"=el"
    exec el.name+"=el"
D = elements.D
T = elements.T
exec D.name+"=D"
exec T.name+"=T"
# Import all elements on "from core import *"
__all__ = ([el.symbol for el in elements]
           + [el.name for el in elements]
           + [D.symbol, D.name, T.symbol, T.name]
           + ['elements'])
