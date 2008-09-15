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

    import danse.gammatable  # loads gammaray data
    from elements import *   # loads symbols for periodic_table, H, He, ...

    for el in periodic_table: print el.symbol, el.gammadata

To implement this, you will need the following in gammatable.py:

    from elements import periodic_table
    
    def _init():
        if 'gammadata' in periodic_table.properties: return
        periodic_table.properties.append('gammadata')

        # Set the default, if any
        elements.elements.Element.gamma = None

        # Set the units
        elements.elements.Element.gamma_units = "Ev"

        # Load the data
        for s,data in gamma_table.iteritems():
            el = periodic_table.symbol(s)
            el.gammadata = data

    # Define the data
    gamma_table = dict(
        Si="Silicon gamma values", 
        O="Oxygen gamma values",
        )

    _init()

You can use similar tricks for isotope specific data:

    from elements import periodic_table

    def _init():
        if 'shells' in periodic_table.properties: return
        periodic_table.properties.append('shells')

        # Set the default.  This is required, even if it is only
        # setting it to None.  If the attribute is missing then the
        # isotope data reverts to the element to supply the value,
        # which is almost certainly not what you want.
        elements.elements.Isotope.shells = None

        # Load the data
        for symbol,data in shell_table.iteritems():
            el = periodic_table.symbol(symbol)
            for iso,isodata in data.iteritems():
                el[iso].shells = isodata

    # Define the data
    shell_table = dict(
        Fe={56: "56-Fe shell info",
            58: "58-Fe shell info",
            }
        )

    _init()

Since your data table is likely to be in its own package, you can
put the following in the package/__init__.py:

     import elements
     table = elements.table  # Allow user to say e.g., gamma.table.Si
     
     import gammatable

You can also define attributes on import that are not loaded directly.
For example, if you don't want to load all the isotope information for
shells immediately, then you can use delayed_load in __init__.py:


     # Delayed loading of shell info
     def _load_shell():
         '''
         Electron shell information for isotopes.
     
         T. Student, Tables of Shell Information
         '''
         import shelltable
     elements.elements.delayed_load(['shells'],_load_shell)

The first argument to delayed_load is the list of all attributes that will
be defined when the module is loaded.  The second argument is the loading
function, whose docstring will appear as the attribute description for
each attribute in the first list.
"""
import copy

# Note: importing self so that we can run as __main__; if we don't then
# elements.Element does not match the implicit __main__.Element
import elements

# Note: __all__ will include all the elements and periodic_table; it is
# defined below.

def delayed_load(all_props,loader,element=True,isotope=False):
    """
    Delayed loading of an element property table.  When any of props
    is first accessed the loader will be called to load the associated
    data.  The help string starts out as the help string for the loader
    function.  If it is an isotope property, be sure to set the
    keyword isotope=True.
    """
    def getter(propname):
        """
        Factory which returns a getter for prop
        """
        def getfn(el):
            for p in all_props:
                delattr(elements.Element,p)
            loader()
            return getattr(el,propname)
        return getfn
    if element:
        for p in all_props:
            setattr(elements.Element,p,property(getter(p),doc=loader.__doc__))
    if isotope:
        for p in all_props:
            setattr(elements.Isotope,p,property(getter(p),doc=loader.__doc__))


# Define the element names from the element table.
class _PeriodicTable(object):
    """
    Defines the period table of the elements.

    Individidual elements are accessed either by table[Z], table.Xx,
    table.name('Name'), table.symbol('Xx') or table.isotope('Xx').
    Individual isotopes are addressable by element[mass_number] or
    table.isotope('n-Xx').
    
    For example, the following all retrieve iron:
        table[26]
        table.Fe
        table.symbol('Fe')
        table.name('iron')
        table.isotope('Fe')

    To get iron-56, use:
        table[26][56]
        table.Fe[56]
        table.isotope('56-Fe')
        
    Deuterium and tritium are defined as table.D and table.T.  Some
    neutron properties are available in table[0]

    To show all the elements in the table, use the iterator:
    
        for element in periodic_table:  # lists the element symbols
            print el.symbol,el.number

    Elements have name, number, symbol and isotopes.
    
    Properties can be added to the elements as needed, including mass, nuclear 
    and X-ray scattering cross sections.  See the help(elements) for details.
    """
    def __init__(self):
        """
        Create one element for each entry in the periodic table.
        """
        self.properties = []
        self._element = {}
        for Z,name in element_symbols.iteritems():
            element = elements.Element(name[0].lower(),name[1],Z)
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
            if isinstance(value,(elements.Element,elements.Isotope)):
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
            if isinstance(attr,elements.Element):
                # If no isotope, return the element
                if isotope == 0:
                    return attr
                # If isotope, check that it is valid
                if isotope in attr.isotopes:
                    return attr[isotope]
            elif isinstance(attr,elements.Isotope):
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
    def __repr__(self):
        # Deuterium and Tritium are special
        if 'symbol' in self.__dict__:
            return self.symbol
        else:
            return "%d-%s"%(self.isotope,self.element.symbol)
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
            self._isotopes[number] = elements.Isotope(self,number)
        return self._isotopes[number]
    def __getitem__(self,number):
        return self._isotopes[number]

    def __iter__(self):
        """
        Process the elements in Z order
        """
        for isotope in self.isotopes:
            yield self._isotopes[isotope]

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
periodic_table = _PeriodicTable()
for el in periodic_table:
    exec el.symbol+"=el"
    exec el.name+"=el"
D = periodic_table.D
T = periodic_table.T
exec D.name+"=D"
exec T.name+"=D"
# Import all elements on "from elements import *"
__all__ = ([el.symbol for el in periodic_table]
           + [el.name for el in periodic_table]
           + [D.symbol, D.name, T.symbol, T.name]
           + ['periodic_table'])

def test():
    # Check that we can access element properties
    assert H.name == "hydrogen"
    assert H.symbol == "H"
    assert H.number == 1
    assert helium.symbol == 'He'

    # Check that isotopes work and produce the correct strings and symbols
    O.add_isotope(18)
    assert H[2].symbol == 'D'
    assert H[3].symbol == 'T'
    assert O[18].symbol == 'O'
    assert str(H[2])=='D'
    assert str(H[3])=='T'
    assert str(O[18])=='18-O'

    # Check the for el in elements works and for iso in el works
    elements = tuple(el for el in periodic_table)
    assert elements[0].number == 0
    assert elements[1].number == 1
    isotopes = tuple(iso for iso in O)
    assert isotopes[0].isotope == 18
        
    # Check that table lookup works and fails appropriately
    Fe.add_isotope(56)
    assert periodic_table.symbol('Fe') == Fe
    assert periodic_table.name('iron') == Fe
    assert periodic_table.isotope('Fe') == Fe
    assert periodic_table.isotope('56-Fe') == Fe[56]
    try:
        periodic_table.symbol('Qu')
    except ValueError,msg:
        assert str(msg) == "unknown element Qu"
    try:
        periodic_table.name('Qu')
    except ValueError,msg:
        assert str(msg) == "unknown element Qu"
    try:
        periodic_table.isotope('Qu')
    except ValueError,msg:
        assert str(msg) == "unknown element Qu"

if __name__ == "__main__": test()
